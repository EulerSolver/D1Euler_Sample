#!/bin/bash
git add .
git commit -m "$(TZ=Asia/Shanghai date '+%Y-%m-%d %H:%M:%S Beijing Time')"
git push origin main