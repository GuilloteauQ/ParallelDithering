PROG = make_bin_folder parallel_dithering

all : $(PROG)

CC        =  gcc
CFLAGS    =  -g -Wall -O2
CPPFLAGS  =  -DDEBUG
LDFLAGS   =  -g -O2 -fopenmp
SRC_FOLDER = src
TARGET_FOLDER = bin

make_bin_folder:
	mkdir -p $(TARGET_FOLDER) > /dev/null

parallel_dithering: $(TARGET_FOLDER)/parallel_dithering.o $(TARGET_FOLDER)/Util.o
	$(CC) -o $@ $^ $(LDFLAGS)

$(TARGET_FOLDER)/%.o: $(SRC_FOLDER)/%.c
	$(CC) -o $@ -c $< $(CFLAGS)

clean :
	@rm -f *.o $(PROG) *.pgm
	@rm -rf $(TARGET_FOLDER)
