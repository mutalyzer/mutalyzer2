CREATE DATABASE `mutalyzer`

USE `mutalyzer`;

CREATE TABLE `BatchJob` (
  `JobID` char(20) NOT NULL,
  `Filter` char(20) NOT NULL,
  `EMail` char(255) NOT NULL,
  `FromHost` char(255) NOT NULL,
  `JobType` char(20) DEFAULT NULL,
  `Arg1` char(20) DEFAULT NULL,
  PRIMARY KEY (`JobID`)
);

CREATE TABLE `BatchQueue` (
  `QueueID` int(5) NOT NULL AUTO_INCREMENT,
  `JobID` char(20) NOT NULL,
  `Input` char(255) NOT NULL,
  `Flags` char(20) DEFAULT NULL,
  PRIMARY KEY (`QueueID`)
);
CREATE TABLE `GBInfo` (
  `AccNo` char(20) NOT NULL DEFAULT '',
  `GI` char(13) DEFAULT NULL,
  `hash` char(32) NOT NULL DEFAULT '',
  `ChrAccVer` char(20) DEFAULT NULL,
  `ChrStart` int(12) DEFAULT NULL,
  `ChrStop` int(12) DEFAULT NULL,
  `orientation` int(2) DEFAULT NULL,
  `url` char(255) DEFAULT NULL,
  PRIMARY KEY (`AccNo`),
  UNIQUE KEY `hash` (`hash`),
  UNIQUE KEY `alias` (`GI`)
);

CREATE TABLE `Link` (
  `mrnaAcc` char(20) NOT NULL,
  `protAcc` char(20) NOT NULL,
  PRIMARY KEY (`mrnaAcc`),
  UNIQUE KEY `protAcc` (`protAcc`)
);

CREATE TABLE `mm1` (
  `hg18` char(50) DEFAULT NULL,
  `hg19` char(50) DEFAULT NULL
);
CREATE TABLE `mm2` (
  `hg18` char(50) DEFAULT NULL,
  `hg19` char(50) DEFAULT NULL
);
