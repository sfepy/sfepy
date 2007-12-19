#!/bin/sh
# 27.08.2007, c

NAME=`basename $0 | awk --field-separator _ '{ print $2}'`
echo 'killing' $NAME'.py ...'
PID=`ps a | grep $NAME | grep -v $0 | grep -v grep | awk '{print $1}'`

if [ -n "$PID" ]
then
    echo 'pid:' $PID
    kill -9 $PID
else
    echo 'no such process!'
fi
echo 'done.'
