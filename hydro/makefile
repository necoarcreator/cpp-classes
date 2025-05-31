# Configuration
INPUT_DIR = out/csv/test/1
OUTPUT_DIR = out/pics
VIDEO_NAME = simulation1.mp4
MPI_PROCS = 4
FFMPEG_FPS = 4
SCALE_WIDTH = 1920
FFMPEG_PRESET = slow
FFMPEG_CRF = 18

.PHONY: all clean run video help

all: clean run video

clean:
	@echo "Cleaning output directory..."
	@rm -rf $(OUTPUT_DIR)/*.png
	@rm -f $(VIDEO_NAME)
	
run:
	@echo "Running MPI simulation..."
	@mpirun -n $(MPI_PROCS) python3 plot_generator.py $(INPUT_DIR) $(OUTPUT_DIR)

video:
	@echo "Creating video..."
	@ffmpeg -y -framerate $(FFMPEG_FPS) -pattern_type glob -i "$(OUTPUT_DIR)/*.png" \
		-c:v libx264 -preset $(FFMPEG_PRESET) -crf $(FFMPEG_CRF) \
		-pix_fmt yuv420p \
		-vf "scale=$(SCALE_WIDTH):-2:flags=neighbor,pad=ceil(iw/2)*2:ceil(ih/2)*2:color=black" \
		$(VIDEO_NAME)
	@echo "Video created: $(VIDEO_NAME)"

help:
	@echo "Usage:"
	@echo "  make all      - Clean, run simulation and create video"
	@echo "  make clean    - Remove generated plots and video"
	@echo "  make run      - Run MPI simulation"
	@echo "  make video    - Create video from generated plots"
	@echo "  make help     - Show this help"
	@echo ""
	@echo "Configuration:"
	@echo "  INPUT_DIR     = $(INPUT_DIR)"
	@echo "  OUTPUT_DIR    = $(OUTPUT_DIR)"
	@echo "  VIDEO_NAME    = $(VIDEO_NAME)"
	@echo "  MPI_PROCS     = $(MPI_PROCS)"
	@echo "  FFMPEG_FPS    = $(FFMPEG_FPS)"
	@echo "  SCALE_WIDTH   = $(SCALE_WIDTH)"
	@echo "  FFMPEG_PRESET = $(FFMPEG_PRESET)"
	@echo "  FFMPEG_CRF    = $(FFMPEG_CRF)"