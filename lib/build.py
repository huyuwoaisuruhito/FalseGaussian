import os

os.chdir('./lib')
os.system('.\\clean.bat')
for dir in os.listdir('./'):
    if os.path.isdir(dir) and dir[-6:] != '.build':
        print('Enter dir:', dir)
        os.chdir(dir)
        os.system('.\\build.bat')
        os.chdir('..')