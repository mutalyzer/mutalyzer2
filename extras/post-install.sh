#!/bin/bash

# Post-install script for Mutalyzer. Run after the setuptools installation
# (python setup.py install).
#
# Notice: The definitions in this file are quite specific to the standard
# Mutalyzer environment. This consists of a Debian stable (Squeeze) system
# with Apache and Mutalyzer using its mod_wsgi module. Debian conventions are
# used throughout. See the README file for more information.
#
# Usage (from the source root directory):
#   sudo bash extras/post-install.sh
#
# Todo:
# - Copy doc to /usr/share/doc
# - Check if MySQL user/database exist before creating
# - General cleanup

set -e
set -u

COLOR_INFO='\033[32m'
COLOR_WARNING='\033[33m'
COLOR_ERROR='\033[31m'
COLOR_END='\033[0m'

# The 'cd /' is a hack to prevent the mutalyzer package under the current
# directory to be used.
PACKAGE_ROOT=$(cd / && python -c 'import mutalyzer; print mutalyzer.package_root()')
BIN_BATCHD=$(which mutalyzer-batchd)
BIN_CACHE_SYNC=$(which mutalyzer-cache-sync)
BIN_MAPPING_UPDATE=$(which mutalyzer-mapping-update)
BIN_WEBSITE=$(which mutalyzer-website.wsgi)
BIN_SOAP_SERVICE=$(which mutalyzer-soap-service.wsgi)
BIN_JSON_SERVICE=$(which mutalyzer-json-service.wsgi)

if [ ! -e /etc/mutalyzer/config ]; then
    echo -e "${COLOR_INFO}Creating /etc/mutalyzer/config${COLOR_END}"
    mkdir -p /etc/mutalyzer
    cp extras/config.example /etc/mutalyzer/config
    chmod -R u=rwX,go=rX /etc/mutalyzer
else
    echo -e "${COLOR_WARNING}Not touching /etc/mutalyzer/config (it exists)${COLOR_END}"
fi

for USERNAME in $(cut -f 1 -d : /etc/passwd); do
    if [ -d "/home/${USERNAME}" ]; then
        echo -e "${COLOR_INFO}Creating /home/${USERNAME}/.config/mutalyzer/config${COLOR_END}"
        echo -e "${COLOR_INFO}Creating /home/${USERNAME}/.cache/mutalyzer${COLOR_END}"
	# Instead of su, use chown later.
        su $USERNAME -c "mkdir -p /home/$USERNAME/.config/mutalyzer"
        su $USERNAME -c "mkdir -p /home/$USERNAME/.cache/mutalyzer"
        su $USERNAME -c "touch /home/$USERNAME/.config/mutalyzer/config"
        su $USERNAME -c "touch /tmp/mutalyzer-$USERNAME.log"
        cat extras/config.user.example > /home/$USERNAME/.config/mutalyzer/config
        sed -i -e "s@<USERNAME>@${USERNAME}@g" /home/$USERNAME/.config/mutalyzer/config
    fi
done

echo -e "${COLOR_INFO}Touching /var/log/mutalyzer.log${COLOR_END}"
touch /var/log/mutalyzer.log
chown www-data:www-data /var/log/mutalyzer.log
chmod u=rw,go=r /var/log/mutalyzer.log

echo -e "${COLOR_INFO}Touching /var/cache/mutalyzer${COLOR_END}"
mkdir -p /var/cache/mutalyzer
chown -R www-data:www-data /var/cache/mutalyzer
chmod -R u=rwX,go=rX /var/cache/mutalyzer

echo -e "${COLOR_INFO}Creating /etc/init.d/mutalyzer-batchd${COLOR_INFO}"
cp extras/init.d/mutalyzer-batchd /etc/init.d/mutalyzer-batchd
sed -i -e "s@<MUTALYZER_BIN_BATCHD>@${BIN_BATCHD}@g" /etc/init.d/mutalyzer-batchd
chmod u=rwx,go=rx /etc/init.d/mutalyzer-batchd

echo -e "${COLOR_INFO}Installing init scripts${COLOR_END}"
update-rc.d -f mutalyzer-batchd remove
update-rc.d mutalyzer-batchd defaults 98 02

echo -e "${COLOR_INFO}Installing crontab${COLOR_END}"
cp extras/cron.d/mutalyzer-cache-sync /etc/cron.d/mutalyzer-cache-sync
sed -i -e "s@<MUTALYZER_BIN_CACHE_SYNC>@${BIN_CACHE_SYNC}@g" /etc/cron.d/mutalyzer-cache-sync
cp extras/cron.d/mutalyzer-mapping-update /etc/cron.d/mutalyzer-mapping-update
sed -i -e "s@<MUTALYZER_BIN_MAPPING_UPDATE>@${BIN_MAPPING_UPDATE}@g" /etc/cron.d/mutalyzer-mapping-update

echo -e "${COLOR_INFO}Creating /etc/apache2/conf.d/mutalyzer.conf${COLOR_END}"
cp extras/apache/mutalyzer.conf /etc/apache2/conf.d/mutalyzer.conf
sed -i -e "s@<MUTALYZER_BIN_WEBSITE>@${BIN_WEBSITE}@g" -e "s@<MUTALYZER_BIN_SOAP_SERVICE>@${BIN_SOAP_SERVICE}@g" -e "s@<MUTALYZER_BIN_JSON_SERVICE>@${BIN_JSON_SERVICE}@g" -e "s@<MUTALYZER_BIN_BATCHD>@${BIN_BATCHD}@g" /etc/apache2/conf.d/mutalyzer.conf
chmod u=rw,go=r /etc/apache2/conf.d/mutalyzer.conf

echo "You will now be asked for the MySQL root password"

# Create databases
cat << EOF | mysql -u root -p
  CREATE USER mutalyzer;
  CREATE DATABASE mutalyzer;
  CREATE DATABASE hg18;
  CREATE DATABASE hg19;
  GRANT ALL PRIVILEGES ON mutalyzer.* TO mutalyzer;
  GRANT ALL PRIVILEGES ON hg18.* TO mutalyzer;
  GRANT ALL PRIVILEGES ON hg19.* TO mutalyzer;
  FLUSH PRIVILEGES;
EOF

echo -e "${COLOR_INFO}Creating tables in hg18 database${COLOR_END}"

# Create ChrName and Mapping table (hg18)
cat << EOF | mysql -u mutalyzer -D hg18
CREATE TABLE ChrName (
  AccNo char(20) NOT NULL,
  name char(20) NOT NULL,
  PRIMARY KEY (AccNo)
);
CREATE TABLE Mapping (
  gene varchar(255) DEFAULT NULL,
  transcript varchar(20) NOT NULL DEFAULT '',
  version smallint(6) DEFAULT NULL,
  chromosome varchar(40) DEFAULT NULL,
  orientation char(1) DEFAULT NULL,
  start int(11) unsigned DEFAULT NULL,
  stop int(11) unsigned DEFAULT NULL,
  cds_start int(11) unsigned DEFAULT NULL,
  cds_stop int(11) unsigned DEFAULT NULL,
  exon_starts longblob NOT NULL,
  exon_stops longblob NOT NULL,
  protein varchar(20) DEFAULT NULL,
  source varchar(20) DEFAULT NULL,
  INDEX (transcript)
);
INSERT INTO ChrName (AccNo, name) VALUES
('NC_000001.9', 'chr1'),
('NC_000002.10', 'chr2'),
('NC_000003.10', 'chr3'),
('NC_000004.10', 'chr4'),
('NC_000005.8', 'chr5'),
('NC_000006.10', 'chr6'),
('NC_000007.12', 'chr7'),
('NC_000008.9', 'chr8'),
('NC_000009.10', 'chr9'),
('NC_000010.9', 'chr10'),
('NC_000011.8', 'chr11'),
('NC_000012.10', 'chr12'),
('NC_000013.9', 'chr13'),
('NC_000014.7', 'chr14'),
('NC_000015.8', 'chr15'),
('NC_000016.8', 'chr16'),
('NC_000017.9', 'chr17'),
('NC_000018.8', 'chr18'),
('NC_000019.8', 'chr19'),
('NC_000020.9', 'chr20'),
('NC_000021.7', 'chr21'),
('NC_000022.9', 'chr22'),
('NC_000023.9', 'chrX'),
('NC_000024.8', 'chrY'),
('NC_001807.4', 'chrM'),
('NT_113891.1', 'chr6_cox_hap1'),
('NT_113959.1', 'chr22_h2_hap1');
EOF

echo -e "${COLOR_INFO}Populating Mapping table with NCBI data (hg18)${COLOR_END}"

# Populate Mapping table with NCBI data (hg18)
MAPPING=$(mktemp)
wget "ftp://ftp.ncbi.nih.gov/genomes/H_sapiens/ARCHIVE/BUILD.36.3/mapview/seq_gene.md.gz" -O - | zcat > $MAPPING
echo -e "${COLOR_INFO}Importing NCBI mapping data, this may take a few minutes (hg18)${COLOR_END}"
$($BIN_MAPPING_UPDATE hg18 $MAPPING reference)
rm $MAPPING

echo -e "${COLOR_INFO}Creating tables in hg19 database${COLOR_END}"

# Create ChrName and Mapping table (hg19)
cat << EOF | mysql -u mutalyzer -D hg19
CREATE TABLE ChrName (
  AccNo char(20) NOT NULL,
  name char(20) NOT NULL,
  PRIMARY KEY (AccNo)
);
CREATE TABLE Mapping (
  gene varchar(255) DEFAULT NULL,
  transcript varchar(20) NOT NULL DEFAULT '',
  version smallint(6) DEFAULT NULL,
  chromosome varchar(40) DEFAULT NULL,
  orientation char(1) DEFAULT NULL,
  start int(11) unsigned DEFAULT NULL,
  stop int(11) unsigned DEFAULT NULL,
  cds_start int(11) unsigned DEFAULT NULL,
  cds_stop int(11) unsigned DEFAULT NULL,
  exon_starts longblob NOT NULL,
  exon_stops longblob NOT NULL,
  protein varchar(20) DEFAULT NULL,
  source varchar(20) DEFAULT NULL,
  INDEX (transcript)
);
INSERT INTO ChrName (AccNo, name) VALUES
('NC_000001.10', 'chr1'),
('NC_000002.11', 'chr2'),
('NC_000003.11', 'chr3'),
('NC_000004.11', 'chr4'),
('NC_000005.9', 'chr5'),
('NC_000006.11', 'chr6'),
('NC_000007.13', 'chr7'),
('NC_000008.10', 'chr8'),
('NC_000009.11', 'chr9'),
('NC_000010.10', 'chr10'),
('NC_000011.9', 'chr11'),
('NC_000012.11', 'chr12'),
('NC_000013.10', 'chr13'),
('NC_000014.8', 'chr14'),
('NC_000015.9', 'chr15'),
('NC_000016.9', 'chr16'),
('NC_000017.10', 'chr17'),
('NC_000018.9', 'chr18'),
('NC_000019.9', 'chr19'),
('NC_000020.10', 'chr20'),
('NC_000021.8', 'chr21'),
('NC_000022.10', 'chr22'),
('NC_000023.10', 'chrX'),
('NC_000024.9', 'chrY'),
('NC_012920.1', 'chrM'),
('NT_167244.1', 'chr6_apd_hap1'),
('NT_113891.2', 'chr6_cox_hap2'),
('NT_167245.1', 'chr6_dbb_hap3'),
('NT_167246.1', 'chr6_mann_hap4'),
('NT_167247.1', 'chr6_mcf_hap5'),
('NT_167248.1', 'chr6_qbl_hap6'),
('NT_167249.1', 'chr6_ssto_hap7'),
('NT_167250.1', 'chr4_ctg9_hap1'),
('NT_167251.1', 'chr17_ctg5_hap1');
EOF

echo -e "${COLOR_INFO}Populating Mapping table with NCBI data (hg19)${COLOR_END}"

# Populate Mapping table with UCSC data (hg19)
MAPPING=$(mktemp)
#wget "ftp://ftp.ncbi.nih.gov/genomes/H_sapiens/mapview/seq_gene.md.gz" -O - | zcat > $MAPPING
wget "ftp://ftp.ncbi.nih.gov/genomes/H_sapiens/ARCHIVE/BUILD.37.2/mapview/seq_gene.md.gz" -O - | zcat > $MAPPING
echo -e "${COLOR_INFO}Importing NCBI mapping data, this may take a few minutes (hg19)${COLOR_END}"
$($BIN_MAPPING_UPDATE hg19 $MAPPING 'GRCh37.p2-Primary Assembly')
rm $MAPPING

echo -e "${COLOR_INFO}Creating tables in mutalyzer database${COLOR_END}"

# Create mutalyzer tables
cat << EOF | mysql -u mutalyzer -D mutalyzer
CREATE TABLE BatchJob (
  JobID char(20) NOT NULL,
  Filter char(20) NOT NULL,
  EMail char(255) NOT NULL,
  FromHost char(255) NOT NULL,
  JobType char(20) DEFAULT NULL,
  Arg1 char(20) DEFAULT NULL,
  PRIMARY KEY (JobID)
);
CREATE TABLE BatchQueue (
  QueueID int(5) NOT NULL AUTO_INCREMENT,
  JobID char(20) NOT NULL,
  Input char(255) NOT NULL,
  Flags char(20) DEFAULT NULL,
  PRIMARY KEY (QueueID),
  KEY JobQueue (JobID,QueueID)
);
CREATE TABLE GBInfo (
  AccNo char(20) NOT NULL DEFAULT '',
  GI char(13) DEFAULT NULL,
  hash char(32) NOT NULL DEFAULT '',
  ChrAccVer char(20) DEFAULT NULL,
  ChrStart int(12) DEFAULT NULL,
  ChrStop int(12) DEFAULT NULL,
  orientation int(2) DEFAULT NULL,
  url char(255) DEFAULT NULL,
  created TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
  PRIMARY KEY (AccNo),
  UNIQUE KEY hash (hash),
  UNIQUE KEY alias (GI),
  INDEX (created)
);
CREATE TABLE Link (
  mrnaAcc char(20) NOT NULL,
  protAcc char(20) NOT NULL,
  PRIMARY KEY (mrnaAcc),
  UNIQUE KEY protAcc (protAcc)
);
EOF

# The remainder is essentially the same as post-upgrade.sh

if [ ! -e /var/www/mutalyzer ]; then
    mkdir -p /var/www/mutalyzer
fi

if [ -e /var/www/mutalyzer/base ]; then
    echo "Removing /var/www/mutalyzer/base"
    rm /var/www/mutalyzer/base
fi

echo -e "${COLOR_INFO}Symlinking /var/www/mutalyzer/base to $PACKAGE_ROOT/templates/base${COLOR_END}"
ln -s $PACKAGE_ROOT/templates/base /var/www/mutalyzer/base

echo "Restarting Apache"
/etc/init.d/apache2 restart

echo "Starting Mutalyzer batch daemon"
/etc/init.d/mutalyzer-batchd start
