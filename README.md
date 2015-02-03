# NestingBlock

This is program of modeling side-chain interaction of peptides in helix/coil states.

It can be easy to calculate the fraction helix of peptides based on Lifson-Roig theory by using following command:

	$ python3 NestingBlock.py test-LifsonRoig

It also can list the equation of equilibrium constant with special force in the block within helix state by using following command:

	$ python3 NestingBlock.py test-NestingBlock

  =========
 | SNOPSIS |
  =========

 $ python3 NestingBlock.py [-p protocol modified] [-c change wnc value] [-v] protocol_file

  ==================
 | FUNCTION LETTERS |
  ==================

-p protocol modified, --protocol=protocol modified

	The <protocol modified> is an python command written in the protocol file. 
	Use this option can over write the the setting already in the protocol file.

-c change wnc value, --chwnc=change wnc value

	<change wnc value> is format like this A.w=1.44, which means change the w value of Ala to 1.44

-v 
	verbose.

  =========
 | Example |
  =========

Calculating the fraction helix of peptie, and overwrite the w value of Ala to 1.2
 
	$ python3 NestingBlock.py -c A.w=1.2 test-LifsonRoig
 
List the equation of equilibrium constant with special force in the block within helix state,
and overwrite the the mode to basic lifson roig (without concering dipole and cap issue)

	$ python3 NestingBlock.py -p mod=\"\" test-NestingBlock
