#!/usr/bin/env python

__doc__ = """

This program requires python 3.6 or higher.

This module, PDBID.py is used to store 82

and 494 PDB IDs. They can be used to test

this program for speed or presence of bugs

"""

__author__     = "Osita Sunday Nnyigide"

__copyright__  = "Copyright 2022, Osita Sunday Nnyigide"

__credits__    = ["Tochukwu Olunna Nnyigide", "Lee Sun-Gu", "Hyun Kyu"]

__license__    = "MIT"

__version__    = "1.2.1"

__maintainer__ = "Osita Sunday Nnyigide"

__email__      = "osita@protein-science.com"

__status__     = "Production"

__date__       = "May 11, 2022"

"""
Below is are two lists, one of 82 PDB IDs
and the other 494 PDB IDs that can be used
to test PRAS. Both are distinct IDs and can
be combined to reach 576 PDB IDs!
"""

_82_pdbs = [

			'1aar', '1afc', '1ah0', '1aj0', '1ajo', '1b9o', '1bd8', '1c5g', '1cec', '1chk', '1cjc', '1cmb',
			'1ctf', '1div', '1dzm', '1ekg', '1ew4', '1f4t', '1fsf', '1gbu', '1go0', '1h09', '1htz', '1hw5',
			'1ic2', '1ipw', '1j2v', '1jji', '1jyd', '1kgb', '1lnw', '1mbe', '1mjc', '1msi', '1npk', '1nr7',
			'1osi', '1oxg', '1p7g', '1pii', '1poh', '1q9e', '1qhe', '1qo5', '1rjb', '1rpo', '1sfp', '1smd',
			'1spp', '1stn',  '1udv', '1uhg','1xyn', '1zug', '2a1d', '2d6c', '2dsp', '2fd7', '2lyz', '2rnt',
			'2sil', '2vh7', '2xbi', '2zta', '2zuq', '3cpa', '3d2a', '3dfq', '3enj', '3jut', '3lel', '3lgf',
			'3ssi', '4ake', '4cha', '4g03', '4gcr', '4oee', '4pep', '5fb6', '5gw9', '5io5',

			]


_494_pdbs = [

			'119l', '153l', '16pk', '19hc', '1a12', '1a28', '1a2p', '1a2z', '1a3a', '1a4i', '1a62', '1a6m', '1a73',
			'1a7s', '1a8d', '1a8e', '1a8i', '1a92', '1aac', '1aay', '1aba', '1ads', '1agj', '1ah7', '1aho', '1aie',
			'1ajj', '1ajs', '1ak0', '1ako', '1akr', '1amf', '1amm', '1amp', '1aoh', '1aop', '1aqb', '1aqu', '1aqz',
			'1arb', '1aru', '1atg', '1atl', '1atz', '1auo', '1avw', '1axn', '1ay7', '1ayl', '1b0t', '1b0y', '1b16',
			'1b3a', '1b4k', '1b4v', '1b5e', '1b67', '1b6a', '1b6g', '1b8o', '1b9w', '1bab', '1bb1', '1bbh', '1bbz',
			'1bd0', '1bdl', '1bdo', '1bec', '1beh', '1ben', '1bf4', '1bf6', '1bfd', '1bfg', '1bg2', '1bg6', '1bgc',
			'1bgf', '1bhp', '1bi5', '1bj7', '1bk0', '1bk7', '1bkb', '1bkj', '1bkr', '1bm8', '1bpi', '1bqc', '1bqk',
			'1brt', '1bs0', '1bs9', '1bsm', '1bte', '1btk', '1bty', '1bu7', '1bu8', '1bue', '1bw9', '1bx4', '1bx7',
			'1bxo', '1byi', '1byq', '1c02', '1c1k', '1c1l', '1c24', '1c3d', '1c3p', '1c3w', '1c52', '1c5e', '1c75',
			'1c90', '1cb0', '1cc8', '1ccz', '1cem', '1ceq', '1cex', '1cf9', '1cgo', '1chd', '1cip', '1cjw', '1cka',
			'1cke', '1cl8', '1cnv', '1cnz', '1cpq', '1cru', '1cs1', '1ctj', '1cv8', '1cvl', '1cxc', '1cxq', '1cxy',
            '1cy4', '1cyd', '1cyo', '1czf', '1czp', '1d2n', '1d3v', '1d7p', '1dbg', '1dbw', '1dci', '1dcs', '1df4',
            '1dfu', '1dgf', '1dhn', '1di6', '1dif', '1din', '1dlf', '1dnl', '1dos', '1doz', '1dp7', '1dps', '1dpt',
			'1dqs', '1dvj', '1dxg', '1eco', '1edg', '1edm', '1egw', '1ek0', '1ek6', '1elk', '1erv', '1erx', '1es5',
			'1etn', '1euw', '1evh', '1ezm', '1fas', '1fdr', '1fds', '1fkj', '1flm', '1flt', '1fmb', '1fna', '1fnc',
			'1ftr', '1fus', '1fxd', '1g3p', '1g82', '1gai', '1gca',	'1gcd', '1gci', '1gd1', '1gof', '1gpe', '1gso',
			'1guq', '1gvp', '1h2r', '1hcl', '1hcr', '1hfc', '1hfb', '1hka', '1hmt', '1hpm', '1htr', '1hxm', '1i4n',
			'1iab', '1ido', '1ifc', '1iib', '1isu', '1ixh', '1jer', '1jet', '1jhg', '1kap', '1koe', '1kp6', '1kpf',
			'1kpt', '1kuh', '1kve', '1l2y', '1lam', '1lbu', '1lcl', '1lkk', '1lst', '1luc', '1m6p', '1mba', '1mct',
			'1mdc', '1mfi', '1mfm', '1mgt', '1mjh', '1mla', '1mml', '1mof', '1mol', '1moq', '1mpg', '1mrj', '1mro',
			'1msk', '1mug', '1mun', '1nar', '1nbc', '1ndd', '1nfn', '1nif', '1nkd', '1nkr', '1nls', '1not', '1nox',
			'1nul', '1nwp', '1nzy', '1oaa', '1onc', '1opd', '1orc', '1osa', '1pcf', '1pda', '1pdo', '1pef', '1pen',
			'1pgs', '1phn', '1plc', '1pmi', '1poa', '1psr', '1ptf', '1pym', '1qau', '1qb7', '1qcx', '1qd1', '1qdd',
			'1qe3', '1qe5', '1qf9', '1qfm', '1qfn', '1qgi', '1qgq', '1qgw', '1qh4', '1qh5', '1qh8', '1qhf', '1qhv',
			'1qip', '1qj4', '1qjd', '1qk5', '1qks', '1ql0', '1qlw', '1qnf', '1qnj', '1qq4', '1qq5', '1qql', '1qqq',
			'1qre', '1qrr', '1qs1', '1qsa', '1qsg', '1qts', '1qtw', '1qu9', '1qup', '1qus', '1ra9', '1rb9', '1rcf',
			'1rge', '1rgg', '1rhs', '1rie', '1rzl', '1sbp', '1slu', '1sml', '1svf', '1svy', '1swu', '1t1d', '1tc1',
			'1tca', '1teq', '1tfe', '1tgs', '1tgx', '1thv', '1tif', '1tml', '1toa', '1tpe', '1tph', '1ttb', '1tud',
			'1tx4', '1tyv', '1uae', '1ubp', '1uch', '1ugi', '1uro', '1ush', '1ute', '1uxy', '1vbo', '1vca', '1vcc',
			'1vfr', '1vfy', '1vhh', '1vie', '1vjs', '1vns', '1vsr', '1w9l', '1wab', '1wap', '1whi', '1xik', '1xjo',
			'1xnb', '1xwl', '1yac', '1yge', '1ytb', '1yve', '1zin', '256b', '2a0b', '2acy', '2ahj', '2arc', '2baa',
			'2bbk', '2bc2', '2bop', '2cba', '2cbp', '2cpg', '2cpl', '2ctc', '2cua', '2dpm', '2dri', '2end', '2eni',
			'2erl', '2fdn', '2gar', '2hbg', '2hft', '2hmz', '2igd', '2ilk', '2izp', '2knt', '2lis', '2mcm', '2mhr',
   			'2msj', '2nac', '2nlr', '2ocj', '2por', '2pth', '2pvb', '2qwc', '2rn2', '2sak', '2sn3', '2spc', '2tgi',
   			'2tnf', '2trc', '2trx', '2xwc', '3bto', '3chb', '3chy', '3cyr', '3ebx', '3eip', '3ezm', '3hts', '3nul',
		    '3oj2', '3pro', '3pte', '3pvi', '3pyp', '3sdh', '3seb', '3sil', '3sti', '3us2', '3vub', '451c', '4cjm',
			'4eug', '4f5s', '4gpd', '4jne', '4lb2', '4lzt', '4n79', '4pga', '4xit', '5cyt', '5hpg', '5nul', '6cel',
			'6gsv', '6ktr', '6m6e', '6n2h', '6o9u', '6v9c', '7a3h', '7atj', '7fd1', '7odc', '7rsa', '8ruc', '9wga',


			]
