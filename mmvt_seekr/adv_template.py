#!/usr/bin/python

'''
A python templating scheme similar to Django's

Usage:  template_object_name = Adv_template('string including $variables to be replaced and commands', {dictionary containing variables and their associated string values})

Like Django's templating system, one can embed commands into the templates, including IF and COMMENT. Commands are surrounded as such: {% command %}.

Examples:

>>> template = "Hi, my name is $NAME and I am $AGE years old."
>>> params = {'NAME':'Johnny', 'AGE':'thirty'}
>>> print(Adv_template(template, params))
Hi, my name is Johnny and I am thirty years old.

>>> template = "{% if $NAME=='Johnny' %}Hi Johnny{% else %}I don't know you{% endif %}"
>>> print(Adv_template(template, params))
Hi Johnny
>>> params = {"NAME":'Susie'}
>>> print(Adv_template(template, params))
I don't know you

'''
import re, string, unittest

class Adv_template():
	def __init__(self, template_string, params): # create a new instance of Advanced template
		self.template_string = template_string
		self.params = params
		self.template_string = self.fix_vars()
		rawoutput = self.parse_commands()
		template_string = string.Template(rawoutput) # create the template
		self.output = template_string.safe_substitute(params) # will substitute every indicated $variable in params. Safely in case template file contains extraneous ($$$) dollar signs

	def get_output(self):
		return self.output

	def __str__(self): # return our output as a string
		return self.output

	def save(self, filename):
		'''saves the template to a file'''
		outfile = open(filename, 'w')
		outfile.writelines(self.output)
		outfile.close()

	def _statement_eval(self, split_command, var_pat):
		# evaluate whether the condition is true
		statement = ' '.join(split_command[1:]) # join the rest of the if statement together to be evaluated
		statement = re.sub(var_pat, r'self.params["\1"]', statement)
		try: # evaluate the condition
			evaluation = eval(statement)
		except KeyError: # if the value has not been set, then assume that it is false
			evaluation = False
		return evaluation

	def parse_commands(self):
		# first, find the blocks & commands
		cmd_pat = re.compile(r'{%(.+?)%}') # the re pattern for matching to a command
		block_pat = re.compile(r'(%\}|^)((\n|.)*?)(\{%|$)') # the re pattern for a block contained by a command
		var_pat = re.compile(r'\$([-_a-zA-Z0-9]+)') # the re pattern to match a variable call in the template (begins with a $)
		commands = re.findall(cmd_pat, self.template_string)
		blocks_tuple = re.findall(block_pat, self.template_string) # returns a 2-d tuple containing all matchings to the 3 regexp groups
		blocks = [n[1] for n in blocks_tuple] # we are only interested in the middle regexp group
		# execute commands of the blocks
		openblock = [] # keeps track of the block that is currently open
		includeblock = [blocks[0]] # always include the first block, no matter what
		including_block = True # we are including the current block
		non_including_depth = 0 # how deep into a non-executed nested block are we?
		skip_block = []
		command_count = 1
		for command in commands: # for every command given
			split_command = command.split() # make command lowercase and split by whitespace
			split_command[0] = split_command[0].lower()
################### the IF command ############################
			if split_command[0] == "if":
				openblock.append(blocks[command_count]) # then we are currently working on this if statement
				skip_block.append(False)
				if non_including_depth == 0: # if we are including the superblock
					evaluation = self._statement_eval(split_command, var_pat)
					if evaluation: # then we can include this block
						including_block = True
						includeblock.append(blocks[command_count]) # include thincludeblock.append(blocks[command_count])is block in the final output
						skip_block[-1] = True
					else: # then we are not including this block
						non_including_depth += 1 # then descend into a depth where we are not going
						including_block = False
				else: # if the superblock is not being included, then increase the non_including_depth
					non_including_depth += 1
					including_block = False

################### ELSEIF command #######################################
			elif (split_command[0] == "elseif") or (split_command[0] == "elif"):
				#if not including_block and non_including_depth == 1 : # if the last statement was false, and we aren't in a non-included superblock
				if not skip_block[-1] and non_including_depth == 1 : # if the last statement was false, and we aren't in a non-included superblock
					openblock.pop()
					openblock.append(blocks[command_count])
					evaluation = self._statement_eval(split_command, var_pat)
					if evaluation: # then we can include this block
						including_block = True
						non_including_depth -= 1
						includeblock.append(blocks[command_count]) # include this block in the final output
						skip_block[-1] = True
					else: # then we are not including this block
						including_block = False
						
				elif including_block:
					non_including_depth += 1
					including_block = False

################### ELSE command #######################################
			elif split_command[0] == "else": # ELSE statement
				#if not including_block and non_including_depth == 1 : # if the last statement was false, and we aren't in a non-included superblock
				if not skip_block[-1] and non_including_depth == 1 : # if the last statement was false, and we aren't in a non-included superblock
					openblock.pop()
					openblock.append(blocks[command_count])
					# then automatically execute this command
					including_block = True
					non_including_depth -= 1
					includeblock.append(blocks[command_count]) # include this block in the final output
				
				elif including_block:
					non_including_depth += 1
					including_block = False
				
			#elif split_command[0] == "for": # FOR loop

################### ENDIF command #######################################
			elif split_command[0] == "endif": # close the IF block
				try:
					openblock.pop()
					skip_block.pop()
				except IndexError:
					print("Unexpected 'endif' command found in namd_input.template.")
					
				if not including_block:
					non_including_depth -= 1 # reduce the depth
					
				assert non_including_depth >= 0, "extra 'endif' statement found!"
				
				if non_including_depth == 0: # then include this block
					including_block = True
					includeblock.append(blocks[command_count])
					#skip_block = []
					
				'''    
				if not including_block: # then we need to see whether ending this block will
					
					if non_including_depth == 0: # then this was the highest superblock that wasn't included
						including_block = True
						includeblock.append(blocks[command_count])
				else:
					includeblock.append(blocks[command_count])
				'''

			#elif split_command[0] == "endfor": # close the FOR loop
################### COMMENT command #######################################
			elif split_command[0] == "comment": 
				openblock.append(blocks[command_count])
				non_including_depth += 1
				including_block = False

################### ENDCOMMENT command #######################################
			elif split_command[0] == "endcomment": # close the comment block
				openblock.pop()
				if non_including_depth > 0: # if we didn't include this block
					non_including_depth -= 1 # reduce the depth
				if non_including_depth == 0: # then this was the highest superblock that wasn't included
					including_block = True
					includeblock.append(blocks[command_count]) # need to include whatever is after this block
			else:
				raise "Unknown Template Command: %s" % (split_command[0],)
			
			command_count += 1
		return ''.join(includeblock)

	def fix_vars(self):
		'''converts all {{ varname }} to $varname in case users prefer that syntax'''
		var_pat = re.compile(r'\{\{\ *(.+)\ *\}\}') # pattern that finds {{ varname }}
		newstring = re.sub(var_pat, r'$\1', self.template_string) # converts to $varname
		return newstring

class File_template(Adv_template):
	'''a template class than replaces values within a string from template file instead of just a string'''
	def __init__(self, template_filename, params):
		self.template_string = ''.join(open(template_filename, 'r').readlines()) # first load the template file and make it a string separated by newlines
		self.params = params
		self.template_string = self.fix_vars()
		rawoutput = self.parse_commands()
		template_string = string.Template(rawoutput) # create the template
		self.output = template_string.safe_substitute(params) # will substitute every indicated $variable in params. Safely in case template file contains extraneous ($$$) dollar signs

	def input_gen(self, filename):
		'''generates input file called (filename) and fills the parameters from the dictionary params'''
		out_file = open(filename, 'w') # open output file
		out_file.write(self.output)
		out_file.close()
