import os
import re
from setuptools import setup, find_packages
from setuptools.command.install import install


__pkg_name__ = 'densecall'
require_file = 'requirements.txt'
package_name = "ch-%s" % __pkg_name__

verstrline = open(os.path.join(__pkg_name__, '__init__.py'), 'r').read()
vsre = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(vsre, verstrline, re.M)
if mo:
    __version__ = mo.group(1)
else:
    raise RuntimeError('Unable to find version string in "{}/__init__.py".'.format(__pkg_name__))


with open(require_file) as f:
    requirements = [r.split()[0] for r in f.read().splitlines()]
try:
    with open('README.md', encoding='utf-8') as f:
        long_description = f.read()
except:
    long_description = ""

setup(
    name=package_name,
    version=__version__,
    packages=['densecall'],
    include_package_data=True,
    install_requires=requirements,
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='linlian',
    author_email='21620151153308@stu.xmu.edu.cn',
    url='https://github.com/n-damo/densecall',
    entry_points = {
        'console_scripts': [
            '{0} = {0}:main'.format(__pkg_name__)
        ]
    },
    
)
