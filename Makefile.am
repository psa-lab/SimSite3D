#
# $Source: /psa/share/repository/pfizer_proj/Makefile.am,v $
# $Revision: 1.3 $
# $Author: vanvoor4 $
# $Date: 2007-02-07 16:15:37 $
#
# $Log: not supported by cvs2svn $
# Revision 1.2  2006/04/18 19:34:57  vanvoor4
# The project is foreign -- ie not GNU style project.
#
# Revision 1.1  2006/03/30 18:05:32  vanvoor4
# Initial checkin
#
#

# Not a GNU style proj
AUTOMAKE_OPTIONS = foreign 

SUBDIRS=src

########### Documentation stuff #####################
EXTRA_DIST = README AUTHORS ChangeLog COPYING INSTALL NEWS

docdir = ${prefix}/share/doc/@PACKAGE@-@VERSION@

install-data-local:
	@$(mkinstalldirs) $(docdir)
	@list='README AUTHORS ChangeLog COPYING INSTALL NEWS'; \
	for p in $$list; do \
	  if test -f $$p; then \
	    echo "$(INSTALL_DATA) $$p $(docdir)/$$p"; \
	    $(INSTALL_DATA) $$p $(docdir)/$$p; \
	else :; fi; \
	done 

install-exec-local:
	@$(mkinstalldirs) $(libdir)/pkgconfig; 
