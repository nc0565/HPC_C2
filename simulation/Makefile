CC=gcc

CFLAGS=-O3 -m64 -std=c99 -Wall -Wextra -Wno-unused-parameter

# -g
# -msse3 -ffast-math # -msse4.2 -mavx2 
# -std=c11
# -march=native

# -cl-finite-math-only -cl-denorms-are-zero

EXE=lbm

PARAM_FILE=../inputs/fast.params

FINAL_STATE_FILE=./final_state.dat
AV_VELS_FILE=./av_vels.dat

REF_FINAL_STATE_FILE=../check/fast.final_state.dat
REF_AV_VELS_FILE=../check/fast.av_vels.dat

DEVICE=0

LDLIBS=-lm

PLATFORM = $(shell uname -s)
ifeq ($(PLATFORM), Darwin)
	LDLIBS += -framework OpenCL
else
	LDLIBS += -lOpenCL
endif

SUBMIT_DIR=submission
SUBMIT_FILES=Makefile lbm.c lbm.h simulation.c utils.c opencl.c kernels.cl
ENV_SCRIPT=env.sh

$(EXE) : utils.o lbm.o simulation.o opencl.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)

%.o: %.c Makefile lbm.h
	$(CC) $(CFLAGS) -o $@ -c $<

# debug: CFLAGS +=  -g
# 	%.o: %.c Makefile lbm.h
# 	$(CC) $(CFLAGS) -o $@ -c $<

run: $(EXE)
	./$(EXE) -a $(AV_VELS_FILE) -f $(FINAL_STATE_FILE) -p $(PARAM_FILE) -d $(DEVICE)

# testVector:
# 	$(CC)  $(CFLAGS) -ftree-vectorizer-verbose=2 -o $@ -c $<

vtune:
	amplxe-cl -quiet -collect hotspots -result-dir r001hs ./$(EXE) -a $(AV_VELS_FILE) -f $(FINAL_STATE_FILE) -p $(PARAM_FILE) -d $(DEVICE)
	amplxe-cl -report hotspots -r r001hs

plot:
	gnuplot final_state.plt

check:
	python ../check/check.py --ref-av-vels-file=$(REF_AV_VELS_FILE) --ref-final-state-file=$(REF_FINAL_STATE_FILE) --av-vels-file=$(AV_VELS_FILE) --final-state-file=$(FINAL_STATE_FILE)

submission:
	@if [ -e $(SUBMIT_DIR) ] ; \
	then \
		echo ; \
		echo "\"$(SUBMIT_DIR)\" already exists and will be deleted" ; \
		read -r -p "Confirm [Y/n] " ; \
		if [[ $$REPLY =~ ^([yY][eE][sS]|[yY]|)$$ ]] ; \
		then \
			rm -rf $(SUBMIT_DIR); \
		else \
			echo Aborting ; \
			false ; \
		fi ; \
	fi

	@rm -rf $(SUBMIT_DIR)
	@mkdir $(SUBMIT_DIR)
	@cp $(SUBMIT_FILES) $(SUBMIT_DIR)
	@if [ -r $(ENV_SCRIPT) ]; then cp $(ENV_SCRIPT) $(SUBMIT_DIR); fi
	@cd $(SUBMIT_DIR) && ../../check/check_submission

.PHONY: clean submission

clean:
	rm -f *.o $(EXE)

