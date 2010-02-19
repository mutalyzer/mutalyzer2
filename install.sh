#!/bin/sh

script="UCSC_update.py"
cron_entry="25 6 \* \* \* python `pwd`/src/$script"

if `echo $0 | grep '/' > /dev/null`; then
  echo "Please run this script from the installation directory."
  exit 1
fi

if ! `crontab -l | grep "$cron_entry" > /dev/null`; then
  echo "Updating cron entry."
  if `crontab -l | grep $script > /dev/null`; then
    echo "Removing old entry."
    crontab -l | grep -v $script | crontab
  fi
  echo "Installing new entry."
  (
    crontab -l
    echo $cron_entry
  ) | crontab
fi
