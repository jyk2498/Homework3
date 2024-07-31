CC = g++

CXXFLAGS = -g
LDFLAGS = -pthread

SRC_DIR = ./src
OBJ_DIR = ./obj

TARGET = RS_code.out

SRCS = $(notdir $(wildcard $(SRC_DIR)/*.cc))
OBJS = $(SRCS:.cc=.o)

OBJECTS = $(patsubst %.o, $(OBJ_DIR)/%.o, $(OBJS))
DEPS = $(OBJECTS:.o=.d)

all: RS_code.out

$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cc
	$(CC) $(CXXFLAGS) -c $< -o $@ -MD $(LDFLAGS)

$(TARGET) : $(OBJECTS)
	$(CC) $(CXXFLAGS) $(OBJECTS) -o $(TARGET) $(LDFLAGS)

.PHONY: clean all

clean:
	rm -f $(OBJECTS) $(DEPS) $(TARGET)

-include $(DEPS)