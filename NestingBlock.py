wncDic = {}
pssDic = {}

def info():
	print('''
 =================
| NestingBlock.py |                 
 =================
 version:1.1
 Date: 2013 Feb 22
 In order to generate the equation by Nesting Block method
 History:
 modified Lifson Roig compelete(LR, Cap)
 dipole concern complete
 add verbose Pss
''')


def usage():
	print('''
 =========	
| SNOPSIS |
 =========
 python NestingBlock.py -p <protocol modified> -c <change wnc value> -[hiv] <protocol_file>

 ==================
| FUNCTION LETTERS |
 ==================

-p <protocol modified>, --protocol=<protocol modified>
	The <protocol modified> is an python command written in the protocol file. 
	Use this option can over write the the setting already in the protocol file.

-c <change wnc value>, --chwnc=<change wnc value>
	<change wnc value> is format like this A.w=1.44, which means change the w value of Ala to 1.44

P.S. All the parameters are CASE SENSITIVE !!!
''')


def StateGenerate( s ):
	length = len(s)
	for i in range(  2 ** length ):
		string = int2bitstr(i)
		string += '0' * ( length - len(string) )
		yield string

def int2bitstr( i ):
	s=''
	while i >1:
		s += str( int(i % 2) )
		i -= i % 2
		i /= 2
	s += str( int( i ) )
	return s

class residue:
	def __init__(self, li):
		self.charge = int(li[0])
		self.w = float(li[1])
		self.n = float(li[2])
		self.c = float(li[3])
		try:
			self.v = float(li[4])
		except:
			self.v = 0.048
			
class pair:
	pair = {}
	def creat( pos1, pos2, var='x'):
		if pos2 < pos1:
			pos1, pos2 = pos2, pos1
		try:
			pair.pair[var]
		except:
			pair.pair[var] = []
		pair.pair[var].append((pos1,pos2))
	
	def check( state ):
		r_dic={}
		for var in pair.pair:
			try:
				r_dic[var]
			except:
				r_dic[var] = 0
			for pos1, pos2 in pair.pair[var]:
				exist = True
				for i in range( pos1, pos2+1):
					if state[i] == '0':
						exist = False
						break
				if exist:
					r_dic[var] += 1
		return r_dic


def loadwnc(filename):
	f = open( filename,"r")
	for line in f:
		line = line.strip()
		if line[0] == '#' or line == "":
			continue
		li = line.strip().split('\t')
		wncDic[ li[0] ] = residue( li[1:])
		
def Hcalc(pattern):
	p = wncDic[pattern[0]].v
	if len(pattern) > 1:
		for i in range(1, len(pattern) - 1):
			p *= wncDic[pattern[i]].w
		p *= wncDic[pattern[-1]].v
	return p

def PssGenerate( sequence ):
	for i in range(len(sequence)):
		for j in range( 1, len(sequence) - i +1):
			try:
				pssDic[sequence[i:i+j]]
			except:
				pssDic[sequence[i:i+j]] = Hcalc(sequence[i:i+j])

def loadprotocol(filename):
	f = open(filename, "r")
	for line in f:
		line = line.strip()
		if line == "" or line[0] =='#':
			continue
		yield( vaidCmd(line) )

def vaidCmd( line ):
	if line.split('=')[0].strip() in ( "sequence", "epsilon", "wncfilename", "mod", "Fhelix_exp", "amidOffSet"):		
		return line
	if "pair.creat" in line:
		return line

def statehelixity( state, denominator ):
	num = 0
	for ch in state:
		if ch == '0':
			continue
		num +=1
	return num / denominator

def LR( sequence, state ):
	''' Lifson-Roig '''
	p = 1
	pos = None
	for i in range(len(sequence)):
		if state[i] == '1' and pos == None:
			pos = i
		elif state[i] == '0' and pos != None:
			p *= pssDic[ sequence[pos:i] ]
			pos = None
	else:
		if pos != None:
			p *= pssDic[ sequence[pos:] ]
	return p
			
def cap( sequence, state ):
	''' modified Lifson-Roig, cap part'''
	if state[0] =='1':
		p = wncDic['ncap'].n
	elif state[1] == '1':
		p = wncDic[sequence[0]].n
	else:
		p = 1
	for i in range(1, len(sequence)-1):
		if state[i]=='0':
			if state[i-1] =='1' and state[i+1] == '1':
				p *= ( wncDic[sequence[i]].n * wncDic[sequence[i]].c ) **0.5
			elif state[i-1] == '1':
				p *= wncDic[sequence[i]].c
			elif state[i+1] == '1':
				p *= wncDic[sequence[i]].n
	if state[-1] == '1':
		p *= wncDic['ccap'].c 
	elif state[-2] == '1':
		p *= wncDic[sequence[-1]].c 
	return p

def poly_gen( var_dic, var):
	i = 0
	while True: 
		if i == 0:
			yield("")
		elif i == 1:
			yield(var)
		else:
			yield(var + str(i) )
		i += 1
		if i > var_dic[var]:
			i =0
			yield( None )

def poly2str(var, ord):
	if ord == 0:
		return("")
	elif ord == 1:
		return var
	else:
		return var + str(ord)
	
class eq:
	def __init__( self, enable=False ): 
		self.enable = enable
		if enable == False:
			self.numerator = 0.0
			self.denominator = 0.0
		else:
			self.numerator = {}
			self.denominator ={}
			var_dic = {} # save the variable and the order
			for var in pair.pair:
				var_dic[var] = len(pair.pair[var])

			var_poly_gen_li =[]
			maxtimes = 1			
			for var in var_dic:
				var_poly_gen_li.append( poly_gen( var_dic, var) )
				maxtimes *= (var_dic[var] + 1)
##			print( [var for var in var_dic])
			field_li = [] 
			for i in range( len(var_poly_gen_li)):
				field_li.append( next(var_poly_gen_li[i]) )
			times =0 
			while times < maxtimes:
				n = 0
				if field_li[0] == None:
					while field_li[n] == None:
						field_li[n] = next( var_poly_gen_li[n] )
						field_li[n+1] = next( var_poly_gen_li[n+1] )
						n += 1
				else:
					key = ''.join(field_li)
					self.numerator[ key ] = 0
					self.denominator[ key ] = 0
					times += 1
					field_li[0] = next( var_poly_gen_li[0] )
			
	def add( self, state, p, h):
		if self.enable == False:
			self.numerator += p * h
			self.denominator += p
		else:
			key =''
			poly = pair.check( state )
			for var in poly:
				key += poly2str( var, poly[var] )
			self.numerator[key] += p * h
			self.denominator[key] += p

	def calc( self, Fhelix_exp=False ):
		if Fhelix_exp == False:
			return self.numerator / self.denominator
		else:
			li = []
			for var in self.numerator:
				li.append( str(self.numerator[var] - self.denominator[var] * Fhelix_exp) + var)
			return ' + '.join(li)

def dipole(epsilon):
	import math
	def prob( epsilon, distance ):
		return math.exp( 3.05877e-8 / ( epsilon * distance ) )
	for seg in pssDic:
		segLen = len(seg)
		for i in range(segLen):
			if wncDic[seg[i]].charge > 0: # Positive residue
				# Repultion with N terminus
				distance = math.sqrt( 2.5e-19 + ( 1.5e-10 * i + 1e-10 ) **2 )
				pssDic[seg] /= prob(epsilon, distance)
				# Attraction with C terminus
				distance = math.sqrt( 2.5e-19 + ( 1.5e-10 * ( segLen - i) +1e-10 ) **2 )
				pssDic[seg] *= prob(epsilon, distance)
			if wncDic[seg[i]].charge < 0: # Negative residue
				# Atrataction with N terminus
				distance = math.sqrt( 2.5e-19 + ( 1.5e-10 * i + 1e-10 ) **2 )
				pssDic[seg] *= prob(epsilon, distance)
				# Repiltion with C terminus
				distance = math.sqrt( 2.5e-19 + ( 1.5e-10 * ( segLen - i) +1e-10 ) **2 )
				pssDic[seg] /= prob(epsilon, distance)
			

def main( sequence, mod="", amidOffSet=-2, Fhelix_exp=False, epsilon=2.0, verbose = False):
	if verbose:
		print( ":: Initialize..")
		if 'c' in mod:
			print("cap concerned (modified lifson-roig)")
		if 'd' in mod:
			print("dipole concerned")
	MaxStateHelixity = len(sequence) + 1 + amidOffSet # number of Amid bond + amidOffSet
	PssGenerate(sequence)
	calc_li =[LR]
	equation = eq(Fhelix_exp)
	if 'c' in mod:
		calc_li.append(cap)
	if 'd' in mod:
		dipole(epsilon) # modify PssDic, concerning dipole

	if verbose:
		print("::Pobability of singles sagament..")
		for key in pssDic:
			print( key, pssDic[key],sep='\t')
		print(":: Generate states..")
		print(sequence, "state probability", sep='\t')
	for state in StateGenerate( sequence ):
		stateprob=1
		for calc in calc_li:
			stateprob *= calc( sequence, state)
		if verbose:
			print ( state, stateprob, sep="\t" )		
		equation.add( state, stateprob, statehelixity( state, MaxStateHelixity ) )
		
	if Fhelix_exp:
		if verbose:
			return equation.calc(Fhelix_exp) + " = 0"
	return equation.calc(Fhelix_exp)


if __name__ == "__main__":
	import sys, getopt
	sequence=""
	mod = "c"
	amidOffSet= -2
	wncfilename="./wnc"
	verbose = False
	epsilon = 2.0
	Fhelix_exp=False
	try:
		opts, args = getopt.getopt(sys.argv[1:], "c:p:hiv", ["chwnc=","protocol=","help", "info", "verbose"])
	except getopt.GetoptError as err: # print help information and exit:
		print( str(err)) # will print something like "option -a not recognized"
		usage()
		sys.exit(2)

	for cmd in loadprotocol( args[0] ):
		exec(cmd)
	for o, a in opts:
		if o in ("-v", "--verbose"):
			verbose = True

	for o, a in opts:
		if o in ("-v", "--verbose"):
			pass
		elif o in ("-h", "--help"):
			usage()
			sys.exit(0)     
		elif o in ("-i", "--info"):
			info()
			sys.exit(0)     
		elif o in ("-p", "--protocol"):
			exec( vaidCmd( a ) )
			if verbose:
				print(":: Change protocol " + a)
		elif o in ( "-c", "--chwnc" ):
			pass
		else:
			assert False, "unhandled option"
	
	loadwnc( wncfilename )
	
	for o, a in opts:
		if o in ( "-c", "--chwnc" ):
			key, value = a.split("=")
			key, atr = key.split('.')
			exec( "wncDic[key]." + atr + "=" + value ) 
			if verbose:
				print(":: Change wnc wncDic[\'" + key + "\']." + atr + "=" + value )
	print( main( sequence, mod, amidOffSet, Fhelix_exp, epsilon, verbose) )
