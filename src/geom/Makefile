include ../Makefile.inc

libSPBgeom: \
	z-AffineGeometry.o \
	z-BVH.o \
	z-lp_solve.o \
	z-predicates.o \
	z-Shapes.o \
	z-ShapeSet.o
	$(AR) rcs libSPBgeom.a z-*.o

z-AffineGeometry.o: AffineGeometry.c AffineGeometry.h
	$(CC) $(CFLAGS) -c AffineGeometry.c -o z-AffineGeometry.o
z-BVH.o: BVH.c BVH.h
	$(CC) $(CFLAGS) -c BVH.c -o z-BVH.o
z-lp_solve.o: lp_solve.c lp_solve.h
	$(CC) $(CFLAGS) -c lp_solve.c -o z-lp_solve.o
z-predicates.o: predicates.c
	$(CC) $(CFLAGS) -c predicates.c -o z-predicates.o
z-Shapes.o: Shapes.c
	$(CC) $(CFLAGS) -c Shapes.c -o z-Shapes.o
z-ShapeSet.o: ShapeSet.c
	$(CC) $(CFLAGS) -c ShapeSet.c -o z-ShapeSet.o

all:
	libSPBgeom

clean:
	rm -f z-*.o
