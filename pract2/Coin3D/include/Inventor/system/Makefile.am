
PublicHeaders = \
	gl.h \
	gl-headers.h \
	inttypes.h
PrivateHeaders = \
	gl-headers.h.in \
	inttypes.h.in
ObsoleteHeaders =

BUILT_SOURCES = \
	gl-headers.h \
	inttypes.h
DISTCLEANFILES = \
	gl-headers.h \
	inttypes.h \
	stamp-h*

gl-headers.h: $(srcdir)/gl-headers.h.in $(top_builddir)/config.status
	@( cd $(top_builddir); $(MAKE) include/Inventor/system/gl-headers.h )

inttypes.h: $(srcdir)/inttypes.h.in $(top_builddir)/config.status
	@( cd $(top_builddir); $(MAKE) include/Inventor/system/inttypes.h )

##$ BEGIN TEMPLATE Make-Include-Common(system, system)
# **************************************************************************
# This template depends on the following variables having been set up before
# its invokation.
#
# PublicHeaders    - header files that should be installed
# PrivateHeaders   - source files that should be distributed, but not
#                    installed - including source-included source files.
# ObsoletedHeaders - header files that should be removed upon installation

# Everything below is autogenerated from a macro, so don't modify
# within the Makefile.am itself. 19991201 mortene.

EXTRA_DIST = $(PrivateHeaders)

if MAC_FRAMEWORK
libsystemincdir = $(MAC_FRAMEWORK_PREFIX)/$(MAC_FRAMEWORK_NAME).framework/Versions/$(MAC_FRAMEWORK_VERSION)/Headers/system
# hack to change libdir value when installing Mac OS X compiler framework
lib@nop@dir = $(MAC_FRAMEWORK_PREFIX)/$(MAC_FRAMEWORK_NAME).framework/Versions/$(MAC_FRAMEWORK_VERSION)/Libraries
else
libsystemincdir = $(includedir)/Inventor/system
endif
libsysteminc_HEADERS = $(PublicHeaders)

install-libsystemincHEADERS: $(libsysteminc_HEADERS)
	@$(NORMAL_INSTALL)
	@if test x"$(libsysteminc_HEADERS)" = x""; then :; else \
	  echo " $(mkinstalldirs) $(DESTDIR)$(libsystemincdir)"; \
	  $(mkinstalldirs) $(DESTDIR)$(libsystemincdir); \
	fi
	@list='$(libsysteminc_HEADERS)'; for p in $$list; do \
	  if test -f "$$p"; then d=; else d="$(srcdir)/"; fi; \
	  f="`echo $$p | sed -e 's|^.*/||'`"; \
	  src="$$d$$p"; dst="$(DESTDIR)$(libsystemincdir)/$$f"; \
	  if cmp -s "$$src" "$$dst"; then :; else \
	    echo " $(INSTALL_HEADER) $$src $$dst"; \
	    $(INSTALL_HEADER) "$$src" "$$dst" || exit 1; \
	  fi \
	done
	@list='$(ObsoletedHeaders)'; for f in $$list; do \
	  file="$(DESTDIR)$(libsystemincdir)/$$f"; \
	  if test -x "$$file"; then \
	    echo " rm -f $$file"; \
	    rm -f "$$file"; \
	  fi \
	done; \
	:

##$ END TEMPLATE Make-Include-Common
