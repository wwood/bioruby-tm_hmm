# bio-tm_hmm

{<img
src="https://secure.travis-ci.org/wwood/bioruby-tm_hmm.png"
/>}[http://travis-ci.org/#!/wwood/bioruby-tm_hmm]

A bioruby plugin for running the transmembrane domain predictor TMHMM automatically on multiple sequences in a FASTA file and manipulation of the results.

## Installation

```
gem install bio-tm_hmm
```

## Usage

```
bio-tm_hmm my.fasta
```

Where my.fasta is a FASTA file with one or more protein sequences in it. Output will be a description of the transmembrane domains predicted by TMHMM.

Other options include -f for printing out the fasta sequences that have some number of transmembrane domains in them, and ignoring those that don't (converse is -g). For instance, to filter out all sequences that have less than 2 predicted transmembrane domains:

```
bio-tm_hmm -f 2 <my.fasta
```

## Developers

To use the library 

```
require 'bio-tm_hmm'
```

The API doc is online. For more code examples see also the test files in
the source tree.
        
## Project home page

Information on the source tree, documentation, issues and how to contribute, see

  http://github.com/wwood/bioruby-tm_hmm

The BioRuby community is on IRC server: irc.freenode.org, channel: #bioruby.

## Cite

  If you use this software, please cite:

[Organellar proteomics reveals hundreds of novel nuclear proteins in the malaria parasite Plasmodium falciparum](genomebiology.com/2012/13/11/R108)

Sophie C Oehring, Ben J Woodcroft, Suzette Moes, Johanna Wetzel, Olivier Dietz, Andreas Pulfer, Chaitali Dekiwadia, Pascal Maeser, Christian Flueck, Kathrin Witmer, Nicolas MB Brancucci, Igor Niederwieser, Paul Jenoe, Stuart A Ralph and Till S Voss	

Genome Biology 2012, 13:R108 doi:10.1186/gb-2012-13-11-r108

## Biogems.info

This Biogem is published at http://biogems.info/index.html#bio-tm_hmm

## Copyright

Copyright (c) 2012 Ben J Woodcroft. See LICENSE.txt for further details.

