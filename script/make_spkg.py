#!/usr/bin/env python
"""Example:
$ ./script/make_spkg.py ../sfepy-2009.3/ /home/share/software/packages/femhub-0.9.3 0 0
 """
import os, sys, shutil
import sfepy

spkg_skeleton_dir = sys.argv[1]
femhub_dir = sys.argv[2]
make_dist = int(sys.argv[3])
install = int(sys.argv[4])

if make_dist:
    os.system('python setup.py sdist')

version = sfepy.__version__
major = version.split()[0]

unpacked_name = os.path.join('dist', 'sfepy-%s' % version)
release_name = unpacked_name + '.tar.gz'

target_dir = os.path.join(spkg_skeleton_dir, 'src')
spkg_name = os.path.basename(os.path.normpath(spkg_skeleton_dir)) + '.spkg'
femhub_target = os.path.join(femhub_dir, 'spkg/standard', spkg_name)

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

if install:
    package_name = os.path.splitext(spkg_name)[0]
    os.remove(os.path.join(femhub_dir, 'spkg/installed', package_name))
    cmd = ' '.join([os.path.join(femhub_dir, 'spd'), '-i', package_name])
    print cmd
    os.system(cmd)
