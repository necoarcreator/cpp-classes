
all: run

run: main.o krootoy.o crest.o tests.o
	g++ main.o krootoy.o crest.o tests.o -o run
	
main.o: main.cpp
	g++ -c main.cpp
	
krootoy.o: krootoy.cpp tests.o
	g++ -c krootoy.cpp
	
crest.o: crest.cpp krootoy.o tests.o
	g++ -c crest.cpp

tests.o: tests.cpp
	g++ -c tests.cpp

clean:
	rm -rf *.o run
	

NUM_RESULTS = 1
ANALYTICS = 2
KORMAK_MODE = 1
TEST = 1
SOLVER_TYPE = 1
FPS = 5

work:
	ifeq ($(SOLVER_TYPE),3)
		ifeq ($(TEST),1)
		  $(MAKE) work_2_k
		else ifeq ($(TEST),0)
		  $(MAKE) work_2_k
		else
		  $(MAKE) work_1_k
		endif
	else
		ifeq ($(TEST),1)
		  $(MAKE) work_2
		else ifeq ($(TEST),0)
		  $(MAKE) work_2
		else
		  $(MAKE) work_1
		endif
	endif
	
work_2: delete_2 main_2 animate_2

work_1: delete_1 main_1 animate_1

work_2_k: delete_2 main_2_k animate_2

work_1_k: delete_1 main_1_k animate_1

delete_2:
	@echo "$(NUM_RESULTS)" | tr ',' '\n' | python3 delete_em_all.py
	@echo "$(ANALYTICS)" | tr ',' '\n' | python3 delete_em_all.py
	
main_2:
	@echo "$(TEST)\n$(SOLVER_TYPE)\n$(NUM_RESULTS)\n$(ANALYTICS)" | ./run

animate_2:
	@echo "$(NUM_RESULTS)\n$(ANALYTICS)" | python3 animation.py
	ffmpeg -y -r $(FPS) -i "results$(NUM_RESULTS)/plot_%03d.png" animation$(NUM_RESULTS).mp4

delete_1:
	@echo "$(NUM_RESULTS)" | tr ',' '\n' | python3 delete_em_all.py
	
main_1:
	@echo "$(TEST)\n$(SOLVER_TYPE)\n$(NUM_RESULTS)" | ./run
	
animate_1:
	@echo "$(NUM_RESULTS)" | python3 animation_solo.py
	ffmpeg -y -r $(FPS) -i "results$(NUM_RESULTS)/plot_%03d.png" animation$(NUM_RESULTS).mp4
	

main_2_k:
	@echo "$(TEST)\n$(SOLVER_TYPE)\n$(KORMAK_MODE)\n$(NUM_RESULTS)\n$(ANALYTICS)" | ./run

main_1_k:
	@echo "$(TEST)\n$(SOLVER_TYPE)\n$(KORMAK_MODE)\n$(NUM_RESULTS)" | ./run

delete:
	@echo "1" | tr ',' '\n' | python3 delete_em_all.py
	@echo "2" | tr ',' '\n' | python3 delete_em_all.py
	@echo "3" | tr ',' '\n' | python3 delete_em_all.py
	@echo "4" | tr ',' '\n' | python3 delete_em_all.py
	@echo "5" | tr ',' '\n' | python3 delete_em_all.py
	@echo "6" | tr ',' '\n' | python3 delete_em_all.py
	@echo "7" | tr ',' '\n' | python3 delete_em_all.py
	@echo "8" | tr ',' '\n' | python3 delete_em_all.py
	@echo "9" | tr ',' '\n' | python3 delete_em_all.py
	
.PHONY: work


