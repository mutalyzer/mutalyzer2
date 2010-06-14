#!/bin/sh

updateCron() {
  cron_entry="$1 python `pwd`/src/$2.py"

  if ! `crontab -l | grep "$cron_entry" > /dev/null`; then
    echo "Updating cron entry."
    if `crontab -l | grep $2 > /dev/null`; then
      echo "Removing old entry."
      crontab -l | grep -v $2 | crontab
    fi
    echo "Installing new entry."
    (
      crontab -l
      echo $cron_entry
    ) | crontab
  fi
}

if `echo $0 | grep '/' > /dev/null`; then
  echo "Please run this script from the installation directory."
  exit 1
fi

updateCron "25 6 \* \* \*" "UCSC_update" 
updateCron "*/1 \* \* \* \*" "BatchChecker"

cat << EOF > .htaccess
SetHandler mod_python
PythonHandler src/handler
PythonPath "sys.path + ['`pwd`/src']"
PythonDebug On

RewriteEngine on
RewriteRule Variant_info.php Variant_info
EOF

chmod go+rx . src src/Modules templates
chmod go+r .htaccess mutalyzer.conf src/*.py src/Modules/*.py templates/*
chmod go+rw var
