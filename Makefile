# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#                   INSTALLING CEGMA FILES
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#     Copyright (C) 2006 -      Genis PARRA FARRE
#                          
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#  INSTALLDIR is the only variable that should be modified
#    (unless you were not using GNUmake, of course, because in such
#     case you must fix all the GNUmake specific commands... ;^D).
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
### PATHS

       BBIN = /bin
       UBIN = /usr/bin
       LBIN = /usr/local/bin
 #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 #  REMEMBER: INSTALLDIR is the only variable that should be modified
 INSTALLDIR = $(LBIN)
 #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#
### MAKE VARS

.PHONY    : clean cleanall forceall
.DEFAULT  : all
.PRECIOUS : %.sh %.awk %.pl %.pm %.ps %.tgz %.gz \
            %.rc %.gff %.gb %.fa README
.SECONDARY :
.SUFFIXES :

#
### INTERPRETERS

       DATE = $(shell date)
       USER = $(shell whoami)

        AWK = $(firstword $(shell which awk))
       AWKV = $(shell $(AWK) --version | \
                $(AWK) '$$0 ~ /^GNU Awk/ { print $$0; }')

       BASH = $(firstword $(shell which bash))
      BASHV = $(shell $(BASH) --version | \
                $(AWK) '$$0 ~ /^GNU bash/ { print $$0; }')

       PERL = $(firstword $(shell which perl))
      PERLV = $(shell $(PERL) --version | \
                $(AWK) '$$0 ~ /^This/ { \
                           gsub(/^This is perl,[ \t]+/,"",$$0); \
                           print $$0; \
                         }') #'
#
### COMMANDS
# You can update the following paths to your system settings

         CP = $(firstword $(shell which cp   )) 
         CHM = $(firstword $(shell which chmod))
#
### LOCAL PATHS

       WORK = .
        SRC = $(WORK)/src
        BIN = $(WORK)/bin

#
### LOCAL VARS

     BINEXT = .sh .awk .pl .pm

# main scripts
      CEGMA = cegma
  LOCAL_MAP = local_map
 GENOME_MAP = genome_map

# parsing scripts

    PARSING =  parsewise hmm_select completeness

# building geneid parameters

 SELF_TRAIN = geneid-train make_paramfile



    SRCCODE = $(CEGMA)  $(GENOME_MAP) $(LOCAL_MAP) \
          $(PARSING) $(SELF_TRAIN)

    BINCODE = $(addprefix $(BIN)/, $(SRCCODE))


#
### MAKE RULES

all: header main trailer

test :
	@echo "### SORRY !!! Tests were not implemented yet...";

install : header installbin trailer

forceall : header cleanbin main trailer

clean : header cleanbin trailer

cleanall : header cleanbin trailer

cleanbin :
	-@$(RM) -f $(wildcard $(BIN)/*);
	rmdir $(BIN);

main : $(BIN) $(BINCODE)
	@$(CHM) 755 $^;

installbin : $(BINCODE)
	@echo "### COPYING BIN FILES TO $(INSTALLDIR)...";
	@$(CP) -p $^ $(INSTALLDIR)/;

#
### FINISHING CODE

$(BIN) :
	mkdir $(BIN);

$(addprefix $(BIN)/, $(CEGMA)) : $(addprefix $(SRC)/, $(CEGMA).pl)
	@echo "### Finishing PERL script from \"$<\" -> \"$@\"" ;
	@( echo "#!$(PERL) -w"; cat $< ) > $@;

$(addprefix $(BIN)/, $(LOCAL_MAP)) : $(addprefix $(SRC)/, $(LOCAL_MAP).pl)
	@echo "### Finishing PERL script from \"$<\" -> \"$@\"" ;
	@( echo "#!$(PERL) -w"; cat $< ) > $@;

$(addprefix $(BIN)/, $(GENOME_MAP)) : $(addprefix $(SRC)/, $(GENOME_MAP).pl)
	@echo "### Finishing PERL script from \"$<\" -> \"$@\"" ;
	@( echo "#!$(PERL) -w"; cat $< ) > $@;

$(addprefix $(BIN)/, parsewise) : $(addprefix $(SRC)/, parsewise.pl)
	@echo "### Finishing PERL script from \"$<\" -> \"$@\"" ;
	@( echo "#!$(PERL) -w"; cat $< ) > $@;

$(addprefix $(BIN)/, hmm_select) : $(addprefix $(SRC)/, hmm_select.pl)
	@echo "### Finishing PERL script from \"$<\" -> \"$@\"" ;
	@( echo "#!$(PERL) -w"; cat $< ) > $@;

$(addprefix $(BIN)/, completeness) : $(addprefix $(SRC)/, completeness.pl)
	@echo "### Finishing PERL script from \"$<\" -> \"$@\"" ;
	@( echo "#!$(PERL) -w"; cat $< ) > $@;

$(addprefix $(BIN)/, geneid-train) : $(addprefix $(SRC)/, geneid-train.pl)
	@echo "### Finishing PERL script from \"$<\" -> \"$@\"" ;
	@( echo "#!$(PERL)"; cat $< ) > $@;

$(addprefix $(BIN)/, make_paramfile) : $(addprefix $(SRC)/, make_paramfile.pl)
	@echo "### Finishing PERL script from \"$<\" -> \"$@\"" ;
	@( echo "#!$(PERL) -w"; cat $< ) > $@;


isexec : $(BINCODE)
	@$(CHM) 755 $<;

#
### INFO RULES

info : header getinfo trailer

getinfo :
	@echo "### BASH ###  $(BASH)  ###  $(BASHV)";
	@echo "### AWK ###  $(AWK)  ###  $(AWKV)";
	@echo "### PERL ###  $(PERL)  ###  $(PERLV)";

header :
	@echo "###";
	@echo "### RUNNING MAKEFILE";
	@echo "###";
	@echo "### $(DATE) -- $(USER)";
	@echo "###";

trailer :
	@echo "###";
	@echo "### MAKEFILE DONE...";
	@echo "###";
