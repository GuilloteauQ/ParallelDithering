PROG = make_bin_folder dithering_bw duplicate_image

all : $(PROG)

CC        =  mpicc
CFLAGS    =  -g -Wall -O2
CPPFLAGS  =  -DDEBUG
LDFLAGS   =  -g -lm -O2
SRC_FOLDER = src
TARGET_FOLDER = bin

make_bin_folder:
	mkdir -p bin

duplicate_image: $(TARGET_FOLDER)/duplicate_image.o $(TARGET_FOLDER)/Util.o
	gcc -o $@ $^ $(LDFLAGS)

dithering_bw: $(TARGET_FOLDER)/dithering_bw.o $(TARGET_FOLDER)/Util.o
	mpicc -o $@ $^ $(LDFLAGS)

$(TARGET_FOLDER)/%.o: $(SRC_FOLDER)/%.c
	$(CC) -o $@ -c $< $(CFLAGS)

clean :
	@rm -f *.o $(PROG) *.pgm
	@rm -rf $(TARGET_FOLDER)
