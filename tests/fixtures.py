"""
Fixtures for unit tests.

Each fixture is defined by a function which when called sets up the fixture.
The order of calling can be important (e.g., fixtures using the database such
as :func:`hg19` must be called after the :func:`database` fixture).
"""


import os
import shutil

from mutalyzer.config import settings
from mutalyzer.db import Base, session
from mutalyzer.db.models import (Assembly, Chromosome, Reference,
                                 TranscriptMapping, TranscriptProteinLink)


#: Reference definitions for use with the :func:`cache` fixture.
REFERENCES = {
    'AB026906.1':   {'filename':   'AB026906.1.gb.bz2',
                     'checksum':   '29b003d5a71af74dc61a92d2ef5cd5d9',
                     'geninfo_id': '5295993'},
    'AL449423.14':  {'filename':   'AL449423.14.gb.bz2',
                     'checksum':   '00a014242818a3b003b4c077af9e10e0',
                     'geninfo_id': '16944057'},
    'NM_000059.3':  {'filename':   'NM_000059.3.gb.bz2',
                     'checksum':   'f93216b3a596adab279ebd7903258548',
                     'geninfo_id': '119395733'},
    'NM_003002.2':  {'filename':   'NM_003002.2.gb.bz2',
                     'checksum':   '990aa672364937335365609617df3050',
                     'geninfo_id': '222352156'},
    'NG_012772.1':  {'filename':   'NG_012772.1.gb.bz2',
                     'checksum':   '163881f00c9c26516d52a4ddb34f941f',
                     'geninfo_id': '256574794',
                     'links':      [('NM_052818', 'NP_438169'),
                                    ('NM_001079691', 'NP_001073159'),
                                    ('NM_000059', 'NP_000050'),
                                    ('NM_001136571', 'NP_001130043')]},
    'AA010203.1':   {'filename':   'AA010203.1.gb.bz2',
                     'checksum':   '57cee03becb77ce68a225b9c844afb24',
                     'geninfo_id': '1471230'},
    'NM_000088.3':  {'filename':   'NM_000088.3.gb.bz2',
                     'checksum':   '5d1f23e3c1799bdb5586c6786b5d5744',
                     'geninfo_id': '110349771'},
    'NM_000143.3':  {'filename':   'NM_000143.3.gb.bz2',
                     'checksum':   'c91799f40fdc0466bf7702af14cf070a',
                     'geninfo_id': '299758401'},
    'NM_002001.2':  {'filename':   'NM_002001.2.gb.bz2',
                     'checksum':   '7fd5aa4fe864fd5193f224fca8cea70d',
                     'geninfo_id': '31317229'},
    'NG_008939.1':  {'filename':   'NG_008939.1.gb.bz2',
                     'checksum':   '114a03e16ad2f63531d796c2fb0d7039',
                     'geninfo_id': '211938431',
                     'links':      [('NM_000532', 'NP_000523')]},
    'NM_000193.2':  {'filename':   'NM_000193.2.gb.bz2',
                     'checksum':   '86d03e1cf38c1387d90116539ea0678f',
                     'geninfo_id': '21071042'},
    'NP_064445.1':  {'filename':   'NP_064445.1.gb.bz2',
                     'checksum':   '33ea9315882b4a9d8c33018a201be2fa',
                     'geninfo_id': '9910526'},
    'L41870.1':     {'filename':   'L41870.1.gb.bz2',
                     'checksum':   '91b1e539a053f731f95d230a06710897',
                     'geninfo_id': '793994'},
    'NG_009105.1':  {'filename':   'NG_009105.1.gb.bz2',
                     'checksum':   'f2579e6c4a8ead4566e485aad493ef7e',
                     'geninfo_id': '216548283',
                     'links':      [('NM_020061', 'NP_064445')]},
    'AF230870.1':   {'filename':   'AF230870.1.gb.bz2',
                     'checksum':   '9fefa34f40d94910edb5de34a3f98910',
                     'geninfo_id': '7739657'},
    'NG_012337.1':  {'filename':   'NG_012337.1.gb.bz2',
                     'checksum':   'ad712f4f225398d2b11b4f08110c70e6',
                     'geninfo_id': '254039638',
                     'links':      [('NM_018195', 'NP_060665'),
                                    ('NM_001082969', 'NP_001076438'),
                                    ('NM_001082970', 'NP_001076439'),
                                    ('NM_003002', 'NP_002993'),
                                    ('NM_012459', 'NP_036591')]},
    'NM_203473.1':  {'filename':   'NM_203473.1',
                     'checksum':   'ec8fbdeda11ef8ec953e4ed39e9a84e5',
                     'geninfo_id': '45439330'},
    'NM_000132.3':  {'filename':   'NM_000132.3.gb.bz2',
                     'checksum':   '94569bee76d7c8b1168e17df4fe1dcb4',
                     'geninfo_id': '192448441'},
    'LRG_1':        {'filename':   'LRG_1.xml.bz2',
                     'checksum':   '5b8f5a39fcd9e3005688eddffd482746'},
    'DMD':          {'accession':  'UD_139015194859',
                     'filename':   'UD_139015194859.gb.bz2',
                     'checksum':   '2cc769c3f636c722142c0aae12662bd4',
                     'links':      [('NM_000109', 'NP_000100'),
                                    ('NM_004006', 'NP_003997'),
                                    ('NM_004009', 'NP_004000'),
                                    ('NM_004010', 'NP_004001'),
                                    ('NM_004007', None),
                                    ('NM_004011', 'NP_004002'),
                                    ('NM_004012', 'NP_004003'),
                                    ('NM_004023', 'NP_004014'),
                                    ('NM_004020', 'NP_004011'),
                                    ('NM_004022', 'NP_004013'),
                                    ('NM_004021', 'NP_004012'),
                                    ('NM_004013', 'NP_004004'),
                                    ('NM_004014', 'NP_004005'),
                                    ('NM_004018', 'NP_004009'),
                                    ('NM_004017', 'NP_004008'),
                                    ('NM_004016', 'NP_004007'),
                                    ('NM_004015', 'NP_004006'),
                                    ('NM_004019', 'NP_004010')]},
    'DPYD':         {'accession':  'UD_139015208095',
                     'filename':   'UD_139015208095.gb.bz2',
                     'checksum':   'b2b9d402a6e43f80ce1e9bbb72a3c0c6',
                     'links':      [('NR_046590', None),
                                    ('XM_005270562', 'XP_005270619'),
                                    ('NM_000110', 'NP_000101'),
                                    ('XM_005270561', 'XP_005270618'),
                                    ('XM_005270563', 'XP_005270620'),
                                    ('XM_005270564', 'XP_005270621'),
                                    ('NM_001160301', 'NP_001153773')]},
    'MARK1':        {'accession':  'UD_139015213982',
                     'filename':   'UD_139015213982.gb.bz2',
                     'checksum':   '0d63a8fe5beddeb793940f6ae194b985',
                     'links':      [('NM_018650', 'NP_061120'),
                                    ('XM_005273133', None),
                                    ('XM_005273134', 'XP_005273191'),
                                    ('XM_005273135', None),
                                    ('XM_005273136', None)]},
    'A1BG':         {'accession':  'UD_139015218717',
                     'filename':   'UD_139015218717.gb.bz2',
                     'checksum':   'e179de8b248806815394c4f7496ba872',
                     'links':      [('NM_001207009', 'NP_001193938'),
                                    ('NM_198458', 'NP_940860'),
                                    ('XM_005258578', 'XP_005258635'),
                                    ('XM_005258577', 'XP_005258634'),
                                    ('NR_015380', None),
                                    ('NM_130786', 'NP_570602'),
                                    ('XM_005258393', 'XP_005258450')]},
    'chr9_reverse': {'accession':  'UD_139015349377',
                     'filename':   'UD_139015349377.gb.bz2',
                     'checksum':   'd21f92d09116c4831ce8d3ef832aa281',
                     'links':      [('NM_001195250', 'NP_001182179'),
                                    ('NR_036576', None),
                                    ('NR_036577', None),
                                    ('NM_001195252', 'NP_001182181'),
                                    ('NM_001195248', 'NP_001182177'),
                                    ('NM_175069', 'NP_778239'),
                                    ('NM_175073', 'NP_778243'),
                                    ('NM_001195251', 'NP_001182180'),
                                    ('NM_001195254', 'NP_001182183'),
                                    ('NR_036578', None),
                                    ('NR_036579', None),
                                    ('NM_001195249', 'NP_001182178')]},
    'COL1A1':       {'accession':   'UD_139022298843',
                     'filename':    'UD_139022298843.gb.bz2',
                     'checksum':    '815517e36fb380b52842ace6a6e78637',
                     'links':       [('XM_005257059', 'XP_005257116'),
                                     ('XM_005257058', 'XP_005257115'),
                                     ('NM_000088', 'NP_000079')]}}


def database():
    """
    Fixture for database table definitions.
    """
    Base.metadata.create_all(session.get_bind())
    session.commit()


def hg19():
    """
    Fixture for GRCh37/hg19 genome assembly with chromosomes.
    """
    assembly = Assembly('GRCh37', 9606, 'Homo sapiens', alias='hg19')
    session.add(assembly)

    session.add_all(Chromosome(assembly, name, accession, organelle_type)
                    for accession, name, organelle_type in [
            ('NC_000001.10', 'chr1', 'chromosome'),
            ('NC_000002.11', 'chr2', 'chromosome'),
            ('NC_000003.11', 'chr3', 'chromosome'),
            ('NC_000004.11', 'chr4', 'chromosome'),
            ('NC_000005.9', 'chr5', 'chromosome'),
            ('NC_000006.11', 'chr6', 'chromosome'),
            ('NC_000007.13', 'chr7', 'chromosome'),
            ('NC_000008.10', 'chr8', 'chromosome'),
            ('NC_000009.11', 'chr9', 'chromosome'),
            ('NC_000010.10', 'chr10', 'chromosome'),
            ('NC_000011.9', 'chr11', 'chromosome'),
            ('NC_000012.11', 'chr12', 'chromosome'),
            ('NC_000013.10', 'chr13', 'chromosome'),
            ('NC_000014.8', 'chr14', 'chromosome'),
            ('NC_000015.9', 'chr15', 'chromosome'),
            ('NC_000016.9', 'chr16', 'chromosome'),
            ('NC_000017.10', 'chr17', 'chromosome'),
            ('NC_000018.9', 'chr18', 'chromosome'),
            ('NC_000019.9', 'chr19', 'chromosome'),
            ('NC_000020.10', 'chr20', 'chromosome'),
            ('NC_000021.8', 'chr21', 'chromosome'),
            ('NC_000022.10', 'chr22', 'chromosome'),
            ('NC_000023.10', 'chrX', 'chromosome'),
            ('NC_000024.9', 'chrY', 'chromosome'),
            ('NT_167244.1', 'chr6_apd_hap1', 'chromosome'),
            ('NT_113891.2', 'chr6_cox_hap2', 'chromosome'),
            ('NT_167245.1', 'chr6_dbb_hap3', 'chromosome'),
            ('NT_167246.1', 'chr6_mann_hap4', 'chromosome'),
            ('NT_167247.1', 'chr6_mcf_hap5', 'chromosome'),
            ('NT_167248.1', 'chr6_qbl_hap6', 'chromosome'),
            ('NT_167249.1', 'chr6_ssto_hap7', 'chromosome'),
            ('NT_167250.1', 'chr4_ctg9_hap1', 'chromosome'),
            ('NT_167251.1', 'chr17_ctg5_hap1', 'chromosome'),
            ('NC_012920.1', 'chrM', 'mitochondrion')])

    session.commit()


def hg19_transcript_mappings():
    """
    Fixture for some selected transcript mappings in the GRCh37/hg19 genome
    assembly. Depends on the :func:`hg19` fixture.
    """
    chromosome_1 = Chromosome.query.filter_by(accession='NC_000001.10').one()
    chromosome_3 = Chromosome.query.filter_by(accession='NC_000003.11').one()
    chromosome_6 = Chromosome.query.filter_by(accession='NC_000006.11').one()
    chromosome_7 = Chromosome.query.filter_by(accession='NC_000007.13').one()
    chromosome_8 = Chromosome.query.filter_by(accession='NC_000008.10').one()
    chromosome_11 = Chromosome.query.filter_by(accession='NC_000011.9').one()
    chromosome_20 = Chromosome.query.filter_by(accession='NC_000020.10').one()
    chromosome_22 = Chromosome.query.filter_by(accession='NC_000022.10').one()
    chromosome_x = Chromosome.query.filter_by(accession='NC_000023.10').one()
    chromosome_mt = Chromosome.query.filter_by(accession='NC_012920.1').one()

    session.add_all([chromosome_1, chromosome_6, chromosome_8, chromosome_11,
                     chromosome_20, chromosome_22, chromosome_mt])

    session.add(TranscriptMapping(
            chromosome_11,
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
    session.add(TranscriptMapping(
            chromosome_11,
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
    session.add(TranscriptMapping(
            chromosome_6,
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
    session.add(TranscriptMapping(
            chromosome_22,
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
    session.add(TranscriptMapping(
            chromosome_22,
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
    session.add(TranscriptMapping(
            chromosome_1,
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
    session.add(TranscriptMapping(
            chromosome_1,
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
    session.add(TranscriptMapping(
            chromosome_1,
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
    session.add(TranscriptMapping(
            chromosome_20,
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
    session.add(TranscriptMapping(
            chromosome_8,
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
    session.add(TranscriptMapping(
            chromosome_mt,
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
    session.add(TranscriptMapping(
            chromosome_1,
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
    session.add(TranscriptMapping(
            chromosome_7,
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
    session.add(TranscriptMapping(
            chromosome_x,
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
    session.add(TranscriptMapping(
            chromosome_x,
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
    session.add(TranscriptMapping(
            chromosome_x,
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
    session.add(TranscriptMapping(
            chromosome_x,
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
    session.add(TranscriptMapping(
            chromosome_x,
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
    session.add(TranscriptMapping(
            chromosome_3,
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

    session.commit()


def cache(*references):
    """
    Returns a cache fixture for the given references.
    """
    def cache_with_references():
        for reference in references:
            entry = REFERENCES[reference]
            try:
                accession = entry['accession']
            except KeyError:
                accession = reference
            geninfo_id = entry.get('geninfo_id')

            path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                'data',
                                entry['filename'])
            shutil.copy(path, settings.CACHE_DIR)

            session.add(Reference(accession, entry['checksum'],
                                  geninfo_identifier=geninfo_id))

            for transcript, protein in entry.get('links', []):
                session.add(TranscriptProteinLink(transcript, protein))

        session.commit()

    return cache_with_references
