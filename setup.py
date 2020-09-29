from setuptools import setup, find_packages
from retrobiocat_web import __version__

with open('requirements.txt') as f:
  requirements = f.read().splitlines()

setup(
  name = 'retrobiocat_web',
  packages = find_packages(),
  include_package_data=True,
  version = __version__,
  license='',
  description = 'Retrosynthesis',
  author = 'William Finnigan',
  author_email = 'wjafinnigan@gmail.com',
  url = '',
  download_url = '',
  keywords = ['enzyme'],
  install_requires=requirements,
  classifiers=[
    'Development Status :: 3 - Alpha',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 3'],
)