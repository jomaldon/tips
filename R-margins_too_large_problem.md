Some times in R you get this kind of annoying message 

```diff
- Error in plot.new() : figure margins too large
```

This means, as the message said, that your plot exceed the plot device limits.
Plot device is the place where the plot will be located. By default, in RStudio plots are located on the "Plots" panel (bottom right panel) but it could be plotted to another "devices" like a PDF file that you are saving. So, margins limit depends on the "device" limits.

Here you will find some tips to obtain *that* desired graph.

### 1) Adjust the plot windows
Move the limit of the RStudio plot windows to maximize the plotting area.  Just press the maximization button (top right the plot windows) and also move the vertical limit, and finally, before the new plot, execute dev.off() to clear the device cache.

### 2) Adjust the plot margins
Margins are a plot feature thar could be adjusted. On the example, margins will be set to *3*, cutting a part of the white space on borders, but also a part of the text. On this case, a value of *4* will be fine (not shown).

### 3) Saving the plot
When you save the plot, the device is redirected to the new file that is being created. This way, plot limits will be a property of the output device and of the **width** and **height** parameters.

### 4) Using X11() command
There are an alternative "plotting device" on MS Windows R version that could be called with the **x11()** command. This will open a new empty windows where the next plot you execute will be plotted. This plotting device could be resized at your convenience (or just maximized) before the plot to ensure the maximun size.

As you can see, this plotting device give you different ways to export the graph or just to copy & paste to the clipboard.

This feature is out of the box in MS Windows but, following some instructions, could be activated in MAC or Linux.

> MACOS with Homebrew  
> brew cask install xquartz

> LINUX
> 
