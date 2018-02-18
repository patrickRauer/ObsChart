from distutils.core import setup

setup(
    name='ObsCharts',
    version='1.0',
    packages=['ObsChart', 'ObsChart.finding_chart', 'ObsChart.visibility_chart'],
    url='',
    license='GPL',
    author='Patrick Rauer',
    author_email='j.p.rauer@sron.nl',
    description='Charts for observation', requires=['astropy', 'PyAstronomy', 'matplotlib', 'numpy', 'astroquery']
)
