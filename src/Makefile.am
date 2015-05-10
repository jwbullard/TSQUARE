bin_PROGRAMS = tsquare

tsquare_SOURCES = \
	tsquare.h tsquare.cc

tsquare_CPPFLAGS = -I. -I$(top_srcdir)/common -I/usr/include -I/usr/local/include -I/opt/local/include

tsquare_LDADD = $(top_srcdir)/common/libtsquare.a $(CGAL_LDFLAGS) $(BOOST_LDFLAGS)
