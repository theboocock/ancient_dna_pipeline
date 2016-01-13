ifeq ($(OS),Windows_NT)
        CC=g++
        CFLAGS=-I.
        BOOST_PATH ?= C:/boost/include/boost-1_55
        ZLIB_PATH ?= C:/GnuWin32
        DEPS= kseq.h R_codes.h
        CXXFLAGS=-std=gnu++0x
        LDFLAGS=-lm
        LIBS=-lboost_filesystem-mgw46-mt-1_55 -lboost_system-mgw46-mt-1_55 -lboost_regex-mgw46-mt-1_55 -lboost_iostreams-mgw46-mt-1_55 -lpthread -lz
        %.o: %.c $(DEPS)
	        $(CC) -c -o $@ $< $(CFLAGS)

        SolexaQA++: main_win.cpp
	        $(CC) -static main_win.cpp -o SolexaQA++ -I. -I$(BOOST_PATH) -I$(ZLIB_PATH)/include -L$(ZLIB_PATH)/lib -L$(BOOST_PATH)/../../lib $(LDFLAGS) $(CXXFLAGS) $(LIBS)

else
        UNAME_S := $(shell uname -s)
        ifeq ($(UNAME_S),Linux)
                CC=g++
                CFLAGS=-I.
                BOOST_PATH ?= /usr/
                DEPS= kseq.h R_codes.h

                LDFLAGS=-lm
                LIBS=-lz -lboost_filesystem -lboost_system -lboost_regex -lboost_iostreams -lpthread
                %.o: %.c $(DEPS)
	                $(CC) -c -o $@ $< $(CFLAGS)

                SolexaQA++: main.cpp
	            	$(CC) -static main.cpp -o SolexaQA++ -I$(BOOST_PATH)/include $(CFLAGS) -L$(BOOST_PATH)/lib $(LDFLAGS) $(CXXFLAGS) $(LIBS)
        endif
        ifeq ($(UNAME_S),Darwin)
                CC=g++
                DEPS= kseq.h R_codes.h
                BOOST_PATH ?=/usr/local/
                CXXFLAGS= -std=c++11 -stdlib=libc++ -mmacosx-version-min=10.7
                LDFLAGS=-lm
                LIBS=-lz -lpthread

                SolexaQA++: main.cpp
			$(CC) main.cpp -o SolexaQA++ $(BOOST_PATH)lib/libboost_filesystem.a $(BOOST_PATH)lib/libboost_system.a $(BOOST_PATH)lib/libboost_regex.a $(BOOST_PATH)lib/libboost_iostreams.a -I. -I/usr/local/include/boost -I/usr/local/include/ -L/usr/local/lib $(LDFLAGS) $(CXXFLAGS) $(LIBS)

endif
endif
