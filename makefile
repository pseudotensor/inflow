
CFLAGS = -O2
LDFLAGS = -lm
CC= gcc

all: inf_const inf_solv

SRCA = inf_const.c nrutil.c newt.c fmin.c fdjac.c ludcmp.c lubksb.c lnsrch.c
OBJA = inf_const.o nrutil.o newt.o fmin.o fdjac.o ludcmp.o lubksb.o lnsrch.o
inf_const: $(OBJA) $(SRCA) makefile
	$(CC) $(CFLAGS) -o inf_const $(OBJA) $(LDFLAGS)

SRCB = inf_solv.c
OBJB = inf_solv.o
inf_solv: $(OBJB) $(SRCB) makefile
	$(CC) $(CFLAGS) -o inf_solv $(OBJB) $(LDFLAGS)

clean:
	rm *.o inf_const inf_solv t00

