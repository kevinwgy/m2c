LIB = libparser.a
SAR = ar rcs

CCFLAGS = $(CXXFLAGS)

OBJS = Lexer.o \
       Parser.tab.o \
       Dictionary.o \
       AssignerCore.o

$(LIB): $(OBJS)
	$(SAR) $(LIB) $(OBJS)

Lexer.o: Parser.tab.cpp Lexer.C
	$(CXX) $(CXXFLAGS) -c Lexer.C

Parser.tab.o: Parser.tab.cpp
	$(CXX) $(CXXFLAGS) -c Parser.tab.cpp -fpermissive

Parser.tab.cpp: Parser.y
	bison -d -v -p yyCmdf Parser.y
	mv Parser.tab.c Parser.tab.cpp

Lexer.C: Lexer.l
	flex -oLexer.C -PyyCmdf Lexer.l

clean:
	rm -f $(LIB) $(OBJS) Lexer.C Parser.output Parser.tab.cpp Parser.tab.h
