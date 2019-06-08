#!/usr/bin/env python3

import argparse, git, shutil
import tempfile, tarfile, fileinput
import email.utils, datetime, os

parser = argparse.ArgumentParser(description='Build yadedaily packages.')
parser.add_argument("-d", help="target distribution", action="store", dest="dist", default='buster', type=str)
args = parser.parse_args()
dt = datetime.datetime.now()

# Define variables
dirpath = "./../deb_%s/"%args.dist
#dirpath = tempfile.mkdtemp()
dirpathyade = dirpath + '/yadedaily/'

repoups = git.Repo('.')
versiondebian = dt.strftime("%Y%m%d") + "~" + repoups.head.commit.hexsha[0:7] + "~" + args.dist + "1"
tarballname = 'yadedaily_%s.orig.tar.xz'%(versiondebian)

# Create tempdir
os.mkdir(dirpath)
# Copy buildtree into the tmpdir
shutil.copytree('.', dirpathyade, ignore=shutil.ignore_patterns('.git', 'deb'))

# Create tarball
with tarfile.open('%s/%s'%(dirpath,tarballname), mode='w:xz') as out:
    print('Creating tarball... %s'%tarballname)
    out.add(dirpathyade + '/.', arcname='yadedaily', recursive=True)

# Copy debian-directory into the proper place
shutil.copytree('./scripts/ppa_ci/debian', dirpathyade + '/debian/')

with fileinput.FileInput(dirpathyade + '/debian/changelog', inplace=True) as file:
    for line in file:
        print(line.replace("VERSION", versiondebian), end='')
with fileinput.FileInput(dirpathyade + '/debian/changelog', inplace=True) as file:
    for line in file:
        print(line.replace("DISTRIBUTION", args.dist), end='')
with fileinput.FileInput(dirpathyade + '/debian/changelog', inplace=True) as file:
    for line in file:
        print(line.replace("DATE", email.utils.formatdate(localtime=True)), end='')

print (versiondebian)
print (dirpath)
#shutil.rmtree(dirpath)