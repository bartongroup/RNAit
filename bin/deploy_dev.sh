#!/usr/bin/env bash

WEBHOST='dag-dev'
HTMLDIR='/var/www/html/dag.compbio.dundee.ac.uk/'
UWSGIDIR='/var/www/uwsgi/dag.compbio.dundee.ac.uk/RNAit/'

echo "Copying static content"
echo "======================"
echo

rsync -av $RNAIT_ROOT/htdocs/images/* $WEBHOST:$HTMLDIR/images/
rsync -av $RNAIT_ROOT/htdocs/RNAit/* $WEBHOST:$HTMLDIR/RNAit/

#cleancss -02 'all:off;removeDuplicateRules:on' $RNAIT_ROOT/htdocs/css/RNAit.css > /tmp/RNAit.css
#rsync -av /tmp/RNAit.css $WEBHOST:$HTMLDIR/css/
rsync -av $RNAIT_ROOT/htdocs/css/RNAit.css $WEBHOST:$HTMLDIR/css/

rsync -av $RNAIT_ROOT/htdocs/js/RNAit.js $WEBHOST:$HTMLDIR/js/

echo
echo "Setting permissions"
echo "==================="
echo

ssh $WEBHOST "chmod 0755 $HTMLDIR/RNAit"
ssh $WEBHOST "chmod 0744 $HTMLDIR/RNAit/index.html"
ssh $WEBHOST "chmod 0744 $HTMLDIR/images/*"
ssh $WEBHOST "chmod 0744 $HTMLDIR/css/RNAit.css"
ssh $WEBHOST "chmod 0744 $HTMLDIR/js/RNAit.js"

echo
echo "Copying UWSGI application"
echo "========================="
echo

rsync -av $RNAIT_ROOT/databases/* $WEBHOST:$UWSGIDIR/databases/
rsync -av $RNAIT_ROOT/templates/* $WEBHOST:$UWSGIDIR/templates/
rsync -av $RNAIT_ROOT/etc/RNAit.prod.yaml $WEBHOST:$UWSGIDIR/RNAit.yaml
rsync -av $RNAIT_ROOT/uwsgi/RNAit.py $WEBHOST:$UWSGIDIR/
ssh $WEBHOST "touch $UWSGIDIR/reload_RNAit"
ssh $WEBHOST "chmod 0755 $UWSGIDIR"
ssh $WEBHOST "chmod 0755 $UWSGIDIR/templates"
ssh $WEBHOST "chmod 0755 $UWSGIDIR/databases"
ssh $WEBHOST "chmod 0755 $UWSGIDIR/RNAit.py"
ssh $WEBHOST "chmod 0744 $UWSGIDIR/RNAit.yaml"
ssh $WEBHOST "chmod 0744 $UWSGIDIR/templates/*"
ssh $WEBHOST "chmod 0744 $UWSGIDIR/databases/*"
