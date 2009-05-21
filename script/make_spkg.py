#!/usr/bin/env python
import os, sys, shutil
import sfepy

spkg_skeleton_dir = sys.argv[1]
femhub_dir = sys.argv[2]
make_dist = int(sys.argv[3])

if make_dist:
    os.system('python setup.py sdist')

version = sfepy.__version__
major = version.split()[0]

unpacked_name = os.path.join('dist', 'sfepy-%s' % version)
release_name = unpacked_name + '.tar.gz'

target_dir = os.path.join(spkg_skeleton_dir, 'src')
spkg_name = os.path.basename(spkg_skeleton_dir) + '.spkg'
femhub_target = os.path.join(femhub_dir, 'spkg/optional', spkg_name)

print version
print major
print release_name
print target_dir
print spkg_name
print femhub_target

os.system('tar xzf "%s" -C dist' % release_name)

try:
    shutil.rmtree(target_dir)
except:
    pass
os.rename(unpacked_name, target_dir)

os.system('tar cjf "%s" "%s"' % (spkg_name, spkg_skeleton_dir))

try:
    os.remove(femhub_target)
except:
    pass
os.rename(spkg_name, femhub_target)
