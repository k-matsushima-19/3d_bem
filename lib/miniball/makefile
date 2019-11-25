TARGET = libminiball.so

OBJS_CPP = src/Miniball.o src/Miniball_wrapper.o

all: $(TARGET)


$(TARGET): $(OBJS_CPP)
	g++ -fPIC -shared -o $@ $(OBJS_CPP)


%.o: %.cpp
	g++ -fPIC -o $@ -c $<

%.o: %.c
	gcc -fPIC -o $@ -c $<

install:
	@cp $(TARGET) ../../lib/

clean:
	@rm -f $(TARGET) $(OBJS_CPP)
