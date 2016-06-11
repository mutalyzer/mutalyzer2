"""
Test fixtures.
"""


from __future__ import unicode_literals

import os
import shutil

import pytest
import yaml

from mutalyzer.config import settings as _settings
from mutalyzer.output import Output
from mutalyzer.redisclient import client as redis
from mutalyzer.db.models import (Assembly, Chromosome, Reference,
                                 TranscriptMapping)
from mutalyzer import db as _db


@pytest.fixture(autouse=True)
def settings(request, tmpdir):
    cache_dir = unicode(tmpdir.mkdir('cache'))
    log_file = unicode(tmpdir.join('log').ensure())

    redis_uri = request.config.option.redis_uri

    _settings.configure({
        'DEBUG':        False,
        'TESTING':      True,
        'CACHE_DIR':    cache_dir,
        'LOG_FILE':     log_file,
        'DATABASE_URI': None,
        'REDIS_URI':    redis_uri
    })

    if redis_uri is not None:
        redis.flushdb()

    return _settings


@pytest.fixture
def output(settings):
    return Output('test')


@pytest.fixture
def db(request, settings, database_uri):
    settings.configure({
        'DATABASE_URI': database_uri
    })

    # Mutalyzer create tables automatically if we're using an SQLite
    # in-memory database.
    if database_uri != 'sqlite://':
        _db.Base.metadata.drop_all(_db.session.get_bind())
        _db.Base.metadata.create_all(_db.session.get_bind())

    request.addfinalizer(_db.session.remove)

    return _db


@pytest.fixture(scope='session')
def available_references():
    filename = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                            'data', 'references.yml')
    with open(filename) as f:
        return yaml.safe_load(f)


def _add_links(settings, links):
    """
    Add transcript-protein links to the cache.
    """
    for transcript, protein in links:
        if transcript is not None:
            key = 'ncbi:transcript-to-protein:%s' % transcript
            if protein is not None:
                redis.set(key, protein)
                if '.' in transcript:
                    key = key.rsplit('.', 1)[0]
                    redis.set(key, protein.rsplit('.', 1)[0])
            else:
                redis.setex(key,
                            settings.NEGATIVE_LINK_CACHE_EXPIRATION,
                            '')
        if protein is not None:
            key = 'ncbi:protein-to-transcript:%s' % protein
            if transcript is not None:
                redis.set(key, transcript)
                if '.' in protein:
                    key = key.rsplit('.', 1)[0]
                    redis.set(key, transcript.rsplit('.', 1)[0])
            else:
                redis.setex(key,
                            settings.NEGATIVE_LINK_CACHE_EXPIRATION,
                            '')


@pytest.fixture
def references(request, settings, db, available_references):
    try:
        keys = request.param
    except AttributeError:
        return []

    references = []

    for key in keys:
        entry = available_references[key]
        try:
            accession = entry['accession']
        except KeyError:
            accession = key

        # TODO: use pytest basepath or something?
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                            'data',
                            entry['filename'])
        shutil.copy(path, settings.CACHE_DIR)

        references.append(Reference(accession, entry['checksum'], 'upload'))

        _add_links(settings, entry.get('links', []))

    db.session.add_all(references)
    db.session.commit()

    return references


@pytest.fixture
def links(request, settings, db, available_references):
    try:
        links = request.param
    except AttributeError:
        return []

    _add_links(settings, links)
    return links


@pytest.fixture
def hg19(db):
    """
    Fixture for GRCh37/hg19 genome assembly with chromosomes.
    """
    assembly = Assembly('GRCh37', 9606, 'Homo sapiens', alias='hg19')
    db.session.add(assembly)

    db.session.add_all(Chromosome(assembly, name, accession, organelle)
                       for accession, name, organelle in [
        ('NC_000001.10', 'chr1', 'nucleus'),
        ('NC_000002.11', 'chr2', 'nucleus'),
        ('NC_000003.11', 'chr3', 'nucleus'),
        ('NC_000004.11', 'chr4', 'nucleus'),
        ('NC_000005.9', 'chr5', 'nucleus'),
        ('NC_000006.11', 'chr6', 'nucleus'),
        ('NC_000007.13', 'chr7', 'nucleus'),
        ('NC_000008.10', 'chr8', 'nucleus'),
        ('NC_000009.11', 'chr9', 'nucleus'),
        ('NC_000010.10', 'chr10', 'nucleus'),
        ('NC_000011.9', 'chr11', 'nucleus'),
        ('NC_000012.11', 'chr12', 'nucleus'),
        ('NC_000013.10', 'chr13', 'nucleus'),
        ('NC_000014.8', 'chr14', 'nucleus'),
        ('NC_000015.9', 'chr15', 'nucleus'),
        ('NC_000016.9', 'chr16', 'nucleus'),
        ('NC_000017.10', 'chr17', 'nucleus'),
        ('NC_000018.9', 'chr18', 'nucleus'),
        ('NC_000019.9', 'chr19', 'nucleus'),
        ('NC_000020.10', 'chr20', 'nucleus'),
        ('NC_000021.8', 'chr21', 'nucleus'),
        ('NC_000022.10', 'chr22', 'nucleus'),
        ('NC_000023.10', 'chrX', 'nucleus'),
        ('NC_000024.9', 'chrY', 'nucleus'),
        ('NT_167244.1', 'chr6_apd_hap1', 'nucleus'),
        ('NT_113891.2', 'chr6_cox_hap2', 'nucleus'),
        ('NT_167245.1', 'chr6_dbb_hap3', 'nucleus'),
        ('NT_167246.1', 'chr6_mann_hap4', 'nucleus'),
        ('NT_167247.1', 'chr6_mcf_hap5', 'nucleus'),
        ('NT_167248.1', 'chr6_qbl_hap6', 'nucleus'),
        ('NT_167249.1', 'chr6_ssto_hap7', 'nucleus'),
        ('NT_167250.1', 'chr4_ctg9_hap1', 'nucleus'),
        ('NT_167251.1', 'chr17_ctg5_hap1', 'nucleus'),
        ('NC_012920.1', 'chrM', 'mitochondrion')])

    db.session.commit()

    return assembly


@pytest.fixture
def hg19_transcript_mappings(db, hg19):
    """
    Fixture for some selected transcript mappings in the GRCh37/hg19 genome
    assembly.
    """
    db.session.add(TranscriptMapping(
        hg19.chromosomes.filter_by(name='chr11').one(),
        'refseq',
        'NM_003002',
        'SDHD',
        'forward',
        111957571,
        111966518,
        [111957571, 111958581, 111959591, 111965529],
        [111957683, 111958697, 111959735, 111966518],
        'ncbi',
        transcript=1,
        cds=(111957632, 111965694),
        select_transcript=False,
        version=2))
    db.session.add(TranscriptMapping(
        hg19.chromosomes.filter_by(name='chr11').one(),
        'refseq',
        'NM_012459',
        'TIMM8B',
        'reverse',
        111955524,
        111957522,
        [111955524, 111957364],
        [111956186, 111957522],
        'ncbi',
        transcript=1,
        cds=(111956019, 111957492),
        select_transcript=False,
        version=2))
    db.session.add(TranscriptMapping(
        hg19.chromosomes.filter_by(name='chr11').one(),
        'refseq',
        'NR_028383',
        'TIMM8B',
        'reverse',
        111955524,
        111957522,
        [111955524, 111956702, 111957364],
        [111956186, 111957034, 111957522],
        'ncbi',
        transcript=1,
        cds=None,
        select_transcript=False,
        version=1))
    db.session.add(TranscriptMapping(
        hg19.chromosomes.filter_by(name='chr6').one(),
        'refseq',
        'NM_000500',
        'CYP21A2',
        'forward',
        32006082,
        32009419,
        [32006082, 32006499, 32006871, 32007133, 32007323, 32007526,
         32007782, 32008183, 32008445, 32008646],
        [32006401, 32006588, 32007025, 32007234, 32007424, 32007612,
         32007982, 32008361, 32008548, 32009419],
        'ncbi',
        transcript=1,
        cds=(32006200, 32008911),
        select_transcript=False,
        version=5))
    db.session.add(TranscriptMapping(
        hg19.chromosomes.filter_by(name='chr22').one(),
        'refseq',
        'NM_001145134',
        'CPT1B',
        'reverse',
        51007290,
        51017096,
        [51007290, 51007765, 51008005, 51008722, 51009320, 51009587,
         51009804, 51010435, 51010632, 51011304, 51011949, 51012764,
         51012922, 51014464, 51014627, 51015286, 51015753, 51016204,
         51016978],
        [51007510, 51007850, 51008097, 51008835, 51009472, 51009721,
         51009968, 51010551, 51010737, 51011489, 51012144, 51012848,
         51013029, 51014541, 51014764, 51015463, 51015892, 51016363,
         51017096],
        'ncbi',
        transcript=1,
        cds=(51007767, 51016344),
        select_transcript=False,
        version=1))
    db.session.add(TranscriptMapping(
        hg19.chromosomes.filter_by(name='chr22').one(),
        'refseq',
        'NR_021492',
        'LOC100144603',
        'forward',
        51021455,
        51022356,
        [51021455, 51022027],
        [51021752, 51022356],
        'ncbi',
        transcript=1,
        cds=None,
        select_transcript=False,
        version=1))
    db.session.add(TranscriptMapping(
        hg19.chromosomes.filter_by(name='chr1').one(),
        'refseq',
        'NM_001007553',
        'CSDE1',
        'reverse',
        115259538,
        115300624,
        [115259538, 115261234, 115262200, 115263160, 115266504, 115267842,
         115268832, 115269604, 115272879, 115273129, 115275225, 115276353,
         115276610, 115277063, 115279379, 115280092, 115280584, 115282313,
         115292442, 115300546],
        [115260837, 115261366, 115262363, 115263338, 115266623, 115267954,
         115269007, 115269711, 115273043, 115273269, 115275437, 115276478,
         115276738, 115277144, 115279476, 115280184, 115280693, 115282511,
         115292828, 115300624],
        'ncbi',
        transcript=1,
        cds=(115260790, 115282511),
        select_transcript=False,
        version=1))
    db.session.add(TranscriptMapping(
        hg19.chromosomes.filter_by(name='chr1').one(),
        'refseq',
        'NM_001130523',
        'CSDE1',
        'reverse',
        115259538,
        115300671,
        [115259538, 115261234, 115262200, 115263160, 115266504, 115267842,
         115268832, 115269604, 115272879, 115273129, 115275225, 115276353,
         115276610, 115277063, 115279379, 115280584, 115282313, 115284148,
         115292442, 115300546],
        [115260837, 115261366, 115262363, 115263338, 115266623, 115267954,
         115269007, 115269711, 115273043, 115273269, 115275437, 115276478,
         115276738, 115277144, 115279476, 115280693, 115282511, 115284294,
         115292828, 115300671],
        'ncbi',
        transcript=1,
        cds=(115260790, 115284285),
        select_transcript=False,
        version=1))
    db.session.add(TranscriptMapping(
        hg19.chromosomes.filter_by(name='chr1').one(),
        'refseq',
        'NM_002241',
        'KCNJ10',
        'reverse',
        160007257,
        160040051,
        [160007257, 160039812],
        [160012322, 160040051],
        'ncbi',
        transcript=1,
        cds=(160011183, 160012322),
        select_transcript=False,
        version=4))
    db.session.add(TranscriptMapping(
        hg19.chromosomes.filter_by(name='chr20').one(),
        'refseq',
        'NM_001162505',
        'TMEM189',
        'reverse',
        48740274,
        48770335,
        [48740274, 48744512, 48746083, 48747402, 48760039, 48770054],
        [48741716, 48744724, 48746227, 48747484, 48760158, 48770335],
        'ncbi',
        transcript=1,
        cds=(48741595, 48770174),
        select_transcript=False,
        version=1))
    db.session.add(TranscriptMapping(
        hg19.chromosomes.filter_by(name='chr8').one(),
        'refseq',
        'NM_017780',
        'CHD7',
        'forward',
        61591339,
        61779465,
        [61591339, 61653818, 61693559, 61707545, 61712947, 61714087,
         61720776, 61728946, 61732566, 61734349, 61734583, 61735062,
         61736399, 61741222, 61742881, 61748632, 61749376, 61750227,
         61750635, 61754203, 61754406, 61757423, 61757809, 61761074,
         61761610, 61763052, 61763591, 61763821, 61764578, 61765057,
         61765388, 61766922, 61768534, 61769004, 61773463, 61774755,
         61775107, 61777575],
        [61591641, 61655656, 61693989, 61707686, 61713084, 61714152,
         61720831, 61729060, 61732649, 61734486, 61734704, 61735305,
         61736575, 61741365, 61743136, 61748842, 61749571, 61750394,
         61750814, 61754313, 61754611, 61757622, 61757968, 61761163,
         61761713, 61763181, 61763663, 61763878, 61764806, 61765265,
         61766059, 61767082, 61768761, 61769447, 61773684, 61774895,
         61775211, 61779465],
        'ncbi',
        transcript=1,
        cds=(61653992, 61778492),
        select_transcript=False,
        version=2))
    db.session.add(TranscriptMapping(
        hg19.chromosomes.filter_by(name='chrM').one(),
        'refseq',
        'NC_012920',
        'ND4',
        'forward',
        10760,
        12137,
        [10760],
        [12137],
        'reference',
        transcript=1,
        cds=(10760, 12137),
        select_transcript=True,
        version=1))
    db.session.add(TranscriptMapping(
        hg19.chromosomes.filter_by(name='chr1').one(),
        'refseq',
        'NM_002001',
        'FCER1A',
        'forward',
        159259504,
        159278014,
        [159259504, 159272096, 159272644, 159273718, 159275778, 159277538],
        [159259543, 159272209, 159272664, 159273972, 159276035, 159278014],
        'ncbi',
        transcript=1,
        cds=(159272155, 159277722),
        select_transcript=False,
        version=2))
    db.session.add(TranscriptMapping(
        hg19.chromosomes.filter_by(name='chr7').one(),
        'refseq',
        'XM_001715131',
        'LOC100132858',
        'reverse',
        19828,
        36378,
        [19828, 20834, 31060, 32957, 35335, 36224],
        [19895, 21029, 31437, 33107, 35541, 36378],
        'ncbi',
        transcript=1,
        cds=(19828, 36378),
        select_transcript=False,
        version=2))
    db.session.add(TranscriptMapping(
        hg19.chromosomes.filter_by(name='chrX').one(),
        'refseq',
        'NM_004011',
        'DMD',
        'reverse',
        31137345,
        32430371,
        [31137345, 31144759, 31152219, 31164408, 31165392, 31187560,
         31190465, 31191656, 31196049, 31196786, 31198487, 31200855,
         31222078, 31224699, 31227615, 31241164, 31279072, 31341715,
         31366673, 31462598, 31496223, 31497100, 31514905, 31525398,
         31645790, 31676107, 31697492, 31747748, 31792077, 31838092,
         31854835, 31893305, 31947713, 31950197, 31986456, 32235033,
         32305646, 32328199, 32360217, 32361251, 32364060, 32366523,
         32380905, 32382699, 32383137, 32398627, 32404427, 32407618,
         32408188, 32429869, 32430279],
        [31140047, 31144790, 31152311, 31164531, 31165635, 31187718,
         31190530, 31191721, 31196087, 31196922, 31198598, 31201021,
         31222235, 31224784, 31227816, 31241238, 31279133, 31341775,
         31366751, 31462744, 31496491, 31497220, 31515061, 31525570,
         31645979, 31676261, 31697703, 31747865, 31792309, 31838200,
         31854936, 31893490, 31947862, 31950344, 31986631, 32235180,
         32305818, 32328393, 32360399, 32361403, 32364197, 32366645,
         32381075, 32382827, 32383316, 32398797, 32404582, 32407791,
         32408298, 32430030, 32430371],
        'ncbi',
        transcript=1,
        cds=(31140036, 32430326),
        select_transcript=False,
        version=3))
    db.session.add(TranscriptMapping(
        hg19.chromosomes.filter_by(name='chrX').one(),
        'refseq',
        'NM_004019',
        'DMD',
        'reverse',
        31196312,
        31285024,
        [31196312, 31198487, 31200855, 31222078, 31224699, 31227615,
         31241164, 31279072, 31284927],
        [31196922, 31198598, 31201021, 31222235, 31224784, 31227816,
         31241238, 31279133, 31285024],
        'ncbi',
        transcript=1,
        cds=(31196782, 31284946),
        select_transcript=False,
        version=2))
    db.session.add(TranscriptMapping(
        hg19.chromosomes.filter_by(name='chrX').one(),
        'refseq',
        'NM_004007',
        'DMD',
        'reverse',
        31137345,
        33038317,
        [31137345, 31144759, 31152219, 31164408, 31165392, 31187560,
         31190465, 31191656, 31196049, 31196786, 31198487, 31200855,
         31222078, 31224699, 31227615, 31241164, 31279072, 31341715,
         31366673, 31462598, 31496223, 31497100, 31514905, 31525398,
         31645790, 31676107, 31697492, 31747748, 31792077, 31838092,
         31854835, 31893305, 31947713, 31950197, 31986456, 32235033,
         32305646, 32328199, 32360217, 32361251, 32364060, 32366523,
         32380905, 32382699, 32383137, 32398627, 32404427, 32407618,
         32408188, 32429869, 32456358, 32459297, 32466573, 32472779,
         32481556, 32482703, 32486615, 32490281, 32503036, 32509394,
         32519872, 32536125, 32563276, 32583819, 32591647, 32591862,
         32613874, 32632420, 32662249, 32663081, 32715987, 32717229,
         32827610, 32834585, 32841412, 32862900, 32867845, 33038256],
        [31140047, 31144790, 31152311, 31164531, 31165635, 31187718,
         31190530, 31191721, 31196087, 31196922, 31198598, 31201021,
         31222235, 31224784, 31227816, 31241238, 31279133, 31341775,
         31366751, 31462744, 31496491, 31497220, 31515061, 31525570,
         31645979, 31676261, 31697703, 31747865, 31792309, 31838200,
         31854936, 31893490, 31947862, 31950344, 31986631, 32235180,
         32305818, 32328393, 32360399, 32361403, 32364197, 32366645,
         32381075, 32382827, 32383316, 32398797, 32404582, 32407791,
         32408298, 32430030, 32456507, 32459431, 32466755, 32472949,
         32481711, 32482816, 32486827, 32490426, 32503216, 32509635,
         32519959, 32536248, 32563451, 32583998, 32591754, 32591963,
         32613993, 32632570, 32662430, 32663269, 32716115, 32717410,
         32827728, 32834757, 32841504, 32862977, 32867937, 33038317],
        'ncbi',
        transcript=1,
        cds=(31140036, 32834745),
        select_transcript=False,
        version=2))
    db.session.add(TranscriptMapping(
        hg19.chromosomes.filter_by(name='chrX').one(),
        'refseq',
        'NM_203473',
        'PORCN',
        'forward',
        48367371,
        48379202,
        [48367371, 48368172, 48369683, 48370280, 48370714, 48370977,
         48371223, 48372628, 48372913, 48374105, 48374278, 48374449,
         48375571, 48378763],
        [48367491, 48368344, 48369875, 48370323, 48370895, 48371107,
         48371240, 48372753, 48373013, 48374181, 48374341, 48374534,
         48375681, 48379202],
        'ncbi',
        transcript=1,
        cds=(48368209, 48378864),
        select_transcript=False,
        version=1))
    db.session.add(TranscriptMapping(
        hg19.chromosomes.filter_by(name='chrX').one(),
        'refseq',
        'NM_000132',
        'F8',
        'reverse',
        154064063,
        154250998,
        [154064063, 154088707, 154089993, 154091358, 154124352, 154128141,
         154129646, 154130326, 154132181, 154132571, 154133086, 154134695,
         154156846, 154175973, 154182167, 154185232, 154189350, 154194245,
         154194701, 154197606, 154212962, 154215512, 154221211, 154225248,
         154227754, 154250685],
        [154066027, 154088883, 154090141, 154091502, 154124507, 154128226,
         154129717, 154130442, 154132363, 154132799, 154133298, 154134848,
         154159951, 154176182, 154182317, 154185446, 154189443, 154194416,
         154194962, 154197827, 154213078, 154215580, 154221423, 154225370,
         154227875, 154250998],
        'ncbi',
        transcript=1,
        cds=(154065872, 154250827),
        select_transcript=False,
        version=3))
    db.session.add(TranscriptMapping(
        hg19.chromosomes.filter_by(name='chr3').one(),
        'refseq',
        'NM_000249',
        'MLH1',
        'forward',
        37034841,
        37092337,
        [37034841, 37038110, 37042446, 37045892, 37048482, 37050305,
         37053311, 37053502, 37055923, 37058997, 37061801, 37067128,
         37070275, 37081677, 37083759, 37089010, 37090008, 37090395,
         37091977],
        [37035154, 37038200, 37042544, 37045965, 37048554, 37050396,
         37053353, 37053590, 37056035, 37059090, 37061954, 37067498,
         37070423, 37081785, 37083822, 37089174, 37090100, 37090508,
         37092337],
        'ncbi',
        transcript=1,
        cds=(37035039, 37092144),
        select_transcript=False,
        version=3))
    db.session.add(TranscriptMapping(
        hg19.chromosomes.filter_by(name='chr17').one(),
        'lrg',
        'LRG_1',
        'COL1A1',
        'reverse',
        48261457,
        48279000,
        [48261457, 48263139, 48263678, 48264001, 48264376, 48264845, 48265237,
         48265457, 48265891, 48266103, 48266264, 48266529, 48266738, 48267040,
         48267220, 48267362, 48267688, 48267904, 48268178, 48268744, 48269149,
         48269341, 48269836, 48270001, 48270158, 48270355, 48271304, 48271491,
         48271710, 48271934, 48272082, 48272408, 48272593, 48272795, 48272928,
         48273284, 48273516, 48273675, 48273845, 48273978, 48274371, 48274541,
         48275093, 48275310, 48275522, 48275794, 48276587, 48276779, 48276917,
         48277114, 48278772],
        [48263009, 48263381, 48263868, 48264283, 48264483, 48264898, 48265344,
         48265510, 48265998, 48266156, 48266371, 48266636, 48266899, 48267093,
         48267273, 48267469, 48267741, 48267957, 48268285, 48268851, 48269247,
         48269385, 48269889, 48270054, 48270211, 48270408, 48271402, 48271544,
         48271808, 48271987, 48272189, 48272461, 48272691, 48272839, 48273026,
         48273337, 48273560, 48273728, 48273889, 48274031, 48274424, 48274594,
         48275146, 48275363, 48275566, 48275865, 48276688, 48276814, 48276951,
         48277308, 48279000],
        'ebi',
        transcript=1,
        cds=(48262863, 48278874),
        select_transcript=True))
    db.session.add(TranscriptMapping(
        hg19.chromosomes.filter_by(name='chr1').one(),
        'lrg',
        'LRG_348',
        'CR2',
        'forward',
        207627645,
        207663240,
        [207627645, 207639871, 207641872, 207642145, 207642495, 207643040,
         207644085, 207644342, 207644768, 207646117, 207646890, 207647146,
         207647586, 207648169, 207649579, 207651230, 207652602, 207653323,
         207658809, 207662487],
        [207627821, 207640257, 207642060, 207642244, 207642577, 207643447,
         207644261, 207644432, 207644844, 207646524, 207647066, 207647230,
         207647668, 207648561, 207649764, 207651415, 207652625, 207653398,
         207658917, 207663240],
        'ebi',
        transcript=1,
        cds=(207627764, 207658899),
        select_transcript=True))

    db.session.commit()


def with_references(*references):
    """
    Convenience decorator for parameterizing tests with reference fixtures.

    Allows us to write:

        @with_references('NM_004006.1', 'NM_004006.2')
        def test_references():
            pass

    Instead of:

        @pytest.mark.usefixtures('references')
        @pytest.mark.parametrize('references',
                                 [['NM_004006.1', 'NM_004006.2']],
                                 ids=['NM_004006.1,NM_004006.2'],
                                 indirect=True)
        def test_references():
            pass

    """
    def test_with_references(test):
        return pytest.mark.usefixtures('references')(
            pytest.mark.parametrize('references', [references], indirect=True,
                                    ids=[','.join(references)])(test))
    return test_with_references


def with_links(*links):
    """
    Convenience decorator for parameterizing tests with transcript-protein
    link fixtures.

    Allows us to write:

        @with_links(('NM_018650', 'NP_061120'), ('NM_027221', None))
        def test_links():
            pass

    Instead of:

        @pytest.mark.usefixtures('links')
        @pytest.mark.parametrize('links',
                                 [('NM_018650', 'NP_061120'),
                                  ('NM_027221', None)],
                                 ids=['NM_018650/NP_061120,NM_027221/*'],
                                 indirect=True)
        def test_links():
            pass

    """
    def test_with_links(test):
        return pytest.mark.usefixtures('links')(
            pytest.mark.parametrize(
                'links', [links], indirect=True,
                ids=[','.join('/'.join(a or '*' for a in l)
                              for l in links)])(test))
    return test_with_links
