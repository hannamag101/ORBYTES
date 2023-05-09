import setuptools

setuptools.setup(
    name = 'ORBYTES',
    version = '0.1',
    author = 'Hanna Adamski',
    author_email = 'hanna.adamski@yale.edu',
    description = 'ORBYTES provides users the ability to interact with an API interface to extract, visualize, and integrate satellite trajectories actively tracked by space-track.org. In addition to mapping the orbital paths of payloads in regions of LEO to GEO, ORBYTES functions as a platform for users to gain intuition about the beauty and complexity of orbital mechanics.' , 
    packages = ['ORBYTES'],
    python_requires = '>=3',
    install_requires = ['numpy', 'matplotlib','pandas','configparser','streamlit','mpltools', 'rebound','mpld3', 'Pillow']
)