#include <iostream>

#include <Window.h>

int main()
{
	Window window(800, 600);
	if(!window.isOpened())
	{
		return 1;
	}
	while(window.update())
	{

	}
	return 0;
}
