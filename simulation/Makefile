CC=gcc

CFLAGS=-O3 -g

EXE=lbm

PARAM_FILE=../inputs/input_128x128.params

FINAL_STATE_FILE=./final_state.dat
AV_VELS_FILE=./av_vels.dat

REF_FINAL_STATE_FILE=../check/input_128x128.final_state_ref.dat
REF_AV_VELS_FILE=../check/input_128x128.av_vels_ref.dat

LDLIBS=-lm

$(EXE) : utils.o lbm.o simulation.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)

%.o: %.c Makefile lbm.h
	$(CC) $(CFLAGS) -o $@ -c $<

run: $(EXE)
	./$(EXE) -a $(AV_VELS_FILE) -f $(FINAL_STATE_FILE) -p $(PARAM_FILE)

plot:
	gnuplot final_state.plt

check:
	python ../check/check.py --ref-av-vels-file=$(REF_AV_VELS_FILE) --ref-final-state-file=$(REF_FINAL_STATE_FILE) --av-vels-file=$(AV_VELS_FILE) --final-state-file=$(FINAL_STATE_FILE)

.PHONY: clean
clean:
	rm -f *.o $(EXE)

