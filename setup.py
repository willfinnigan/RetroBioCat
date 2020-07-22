from setuptools import setup, find_packages

setup(
  name = 'retrobiocat_web',
  packages = find_packages(),
  include_package_data=True,
  version = '0.1',
  license='',
  description = 'Retrosynthesis',
  author = 'William Finnigan',
  author_email = 'wjafinnigan@gmail.com',
  url = '',
  download_url = '',
  keywords = ['enzyme'],
  install_requires=['crossrefapi', 'pubchempy', 'networkx', 'matplotlib', 'pyyaml'],
  classifiers=[
    'Development Status :: 3 - Alpha',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 3'],
)