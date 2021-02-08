from setuptools import setup
from re import search as re_search
from os import path as os_path


## README
def read_readme(filename):
    with open(filename, 'r', encoding='utf-8') as f:
        return f.read()


## RELEASE
def get_version(filename):
    with open(filename, 'r') as f:
        m = re_search('"(.+)"', f.readline().split('=')[1])
        if m:
            return m.group(1)

    # with open(filename, 'r') as f:
    #     line = f.readline()
    #     while line:
    #         match = re_search("^## (\d\.\d\.\d)$", line)
    #         if match:
    #             return match.group(1)
    #         line = f.readline()


## EXTRAS INFOS
def get_extras(filename):

    extras_infos = {
        'package'      : '',
        'url'          : '',
        'autors'       : '',
        'corr_authors' : '',
        'descr'  : ''
    }

    with open(filename, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('PACKAGE='):
                extras_infos['package'] = line.splitlines()[0].split('=')[1].lower()
            if line.startswith('URL='):
                extras_infos['url'] = line.splitlines()[0].split('=')[1].lower()
            if line.startswith('AUTHORS='):
                extras_infos['authors'] = line.splitlines()[0].split('=')[1].lower()
            if line.startswith('DESCR='):
                extras_infos['descr'] = line.splitlines()[0].split('=')[1].lower()
            if line.startswith('CORR_AUTHOR='):
                extras_infos['corr_authors'] = line.splitlines()[0].split('=')[1].lower()

    return extras_infos


## PACKAGES NAMES
def get_packages(package):
    with open(os_path.join(package, '.env'), 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('MODULES='):
                modules  = line.splitlines()[0].split('=')[1].lower().split(',')
                packages = [package] + [package+'.'+m for m in modules]
    return packages


# get extras infos
extras_infos = get_extras(os_path.join('extras', '.env'))


## SETUP
setup(
    name                          = extras_infos['package'],
    version                       = get_version(os_path.join(extras_infos['package'], '_version.py')),
    author                        = extras_infos['authors'],
    author_email                  = extras_infos['corr_authors'],
    description                   = extras_infos['descr'],
    long_description              = read_readme('README.md'),
    long_description_content_type = 'text/markdown',
    url                           = extras_infos['url'],
    packages                      = get_packages(extras_infos['package']),
    include_package_data          = True,
    test_suite                    = 'pytest',
    license                       = 'MIT',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires               = '>=3.6',
)
