#! /usr/bin/env python

# Taken from PonyGE
# Copyright (c) 2009-2012 Erik Hemberg and James McDermott
# Hereby licensed under the GNU GPL v3.
# http://ponyge.googlecode.com
#
# but improved/changed a bit since then:

# allow eg n_vars=10 kwarg 
# update to Py3
# represent symbol as string s, not tuple (s, "NT") 
# simplify the derivation using iterators

import sys, copy, re, random, math, operator, pprint, itertools, functools
from .exceptions import InvalidIndividualException
import math

class Grammar(object):
    """Context Free Grammar"""

    def __init__(self, file_name, **kwargs):
        self.rules = {}
        self.non_terminals, self.terminals = set(), set()
        self.start_rule = None

        self.read_bnf_file(file_name, **kwargs)
        # TODO check lcm
        self.production_lcm = math.lcm(*(len(self.rules[k]) for k in self.non_terminals))
        # print('grammar self.production_lcm', self.production_lcm)

    def read_bnf_file(self, file_name, **kwargs):
        """Read a grammar file in BNF format"""
        rule_separator = "::="
        # Don't allow space in NTs, and use lookbehind to match "<"
        # and ">" only if not preceded by backslash. Group the whole
        # thing with capturing parentheses so that split() will return
        # all NTs and Ts. TODO does this handle quoted NT symbols?
        non_terminal_pattern = r"((?<!\\)<\S+?(?<!\\)>)"
        # Use lookbehind again to match "|" only if not preceded by
        # backslash. Don't group, so split() will return only the
        # productions, not the separators.
        production_separator = r"(?<!\\)\|"

        # Read the grammar file
        for line in open(file_name, 'r'):
            if not line.startswith("#") and line.strip() != "":
                # Split rules. Everything must be on one line
                if line.find(rule_separator):
                    lhs, productions = line.split(rule_separator, 1) # 1 split
                    lhs = lhs.strip()
                    if not re.search(non_terminal_pattern, lhs):
                        raise ValueError("lhs is not a NT:", lhs)
                    self.non_terminals.add(lhs)
                    if self.start_rule == None:
                        self.start_rule = lhs
                    # Find terminals and non-terminals
                    tmp_productions = []

                    if productions.strip().startswith("GE_RANGE:"):
                        # Special case: for GE_RANGE:nvars, substitute
                        # 0 | 1 | ... | 9 (assuming nvars=10)
                        varname = productions.strip().split(":")[1]
                        assert varname in kwargs
                        n = kwargs[varname]
                        tmp_productions = []
                        for i in range(n):
                            tmp_production = []
                            symbol = str(i)
                            tmp_production.append(symbol)
                            tmp_productions.append(tmp_production)
                            self.terminals.add(symbol)
                    else:
                        # Usual case: iterate through productions on RHS
                        for production in re.split(production_separator, productions):
                            production = production.strip().replace(r"\|", "|")
                            tmp_production = []
                            for symbol in re.split(non_terminal_pattern, production):
                                symbol = symbol.replace(r"\<", "<").replace(r"\>", ">")
                                if len(symbol) == 0:
                                    continue
                                elif re.match(non_terminal_pattern, symbol):
                                    tmp_production.append(symbol)
                                else:
                                    self.terminals.add(symbol)
                                    tmp_production.append(symbol)

                            tmp_productions.append(tmp_production)
                    # Create a rule
                    if not lhs in self.rules:
                        self.rules[lhs] = tmp_productions
                    else:
                        raise ValueError("lhs should be unique", lhs)
                else:
                    raise ValueError("Each rule must be on one line")

    def __str__(self):
        return "%s %s %s %s" % (self.terminals, self.non_terminals,
                                self.rules, self.start_rule)

def derive_string(grammar, genome):
    """Recursively derive a string given a grammar, genome, and start
    symbol. Track the number of codons used. Don't create the
    derivation tree."""
    s = grammar.start_rule
    used_codons = itertools.count()
    try:
        s = _derive_string(grammar, iter(genome), s, used_codons)
        return s, next(used_codons)
    except StopIteration:
        raise InvalidIndividualException
    
def _derive_string(grammar, genome, s, used_codons):
    if s in grammar.terminals:
        return s
    rule = grammar.rules[s]
    if len(rule) > 1:
        codon = next(genome)
        next(used_codons)
        idx = codon % len(rule)
        prod = rule[idx]
        # print "rule", rule, "codon", codon, "prod", prod
    else:
        prod = rule[0]
        # print "rule", rule, "no codon", "prod", prod
    # recurse 
    return "".join([_derive_string(grammar, genome, s, used_codons)
                   for s in prod])
    

