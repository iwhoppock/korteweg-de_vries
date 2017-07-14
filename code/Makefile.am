ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = \
	Korgi 

Korgi_SOURCES = \
	Korgi.c \
	RKfour.c \
	spatial_finite_difference.c \
	linear_algebra.h 


LDADD = liblinalg.a

lib_LIBRARIES = liblinalg.a

liblinalg_a_SOURCES = \
	linear_algebra.h \
	vector.c \
	vector_add.c

