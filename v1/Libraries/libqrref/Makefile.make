OBJS=zgeqr2.o zgeqrf.o zlarf.o zlarfb.o zlarfg.o zlarft.o zung2r.o zungqr.o  zlapmr.o

all: $(OBJS)
	ar -r libqrref.a $(OBJS)

%.o: %.f
	$(ALF_FC) -c $(ALF_FLAGS_QRREF) $< 

clean: 
	rm -f $(OBJS) libqrref.a
