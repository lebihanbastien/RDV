#include <stdio.h>
#include <stdlib.h>
#include <string.h> // for strchr()


void emptyBuffer()
{
    int c = 0;
    while (c != '\n' && c != EOF)
    {
        c = getchar();
    }
}



int reads(char *input, int length)
{
    char *sp = NULL;
    if (fgets(input, length, stdin) != NULL)
    {
        sp = strchr(input, '\n');

        if (sp != NULL)
        {
            *sp = '\0';
        }
        else
        {
            emptyBuffer();
        }
        return 1;
    }
    else

    {
        emptyBuffer();

        return 0;
    }
}
