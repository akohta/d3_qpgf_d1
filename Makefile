CC      =gcc
CFLAGS  =-O3
LDFLAGS =-lm
SRCDIR  =src
OBJDIR  =src
TARGET1 =example.out
TRGSRC1 =example.c

SRCS=$(wildcard $(SRCDIR)/*.c)
OBJS=$(addprefix $(OBJDIR)/,$(patsubst %.c,%.o,$(notdir $(SRCS)) ))
HEAD=$(wildcard $(SRCDIR)/*.h)

TRGOBJ1=$(OBJS) 

all : $(TARGET1)

$(TARGET1) : $(TRGOBJ1) 
	$(CC) $(LDFLAGS) -o $@ $^ 

$(OBJDIR)/%.o : $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -I $(SRCDIR) -c $< -o $@
        
clean:
	@rm -rf $(TARGET1) $(OBJDIR)/*.o
        
$(OBJS) : $(HEAD)


