# CS-4350 Project

## Bassam Kamal - 367032 | Cityscape - 12/10/2024

![s4.png](https://github.com/Bbrizly/Cityscape/blob/main/data/screenshots/s4.png)


## Controls:
WASD            for movement (forwards, backwards, left, right)
Left Shift      to speed up
Space           for up
Left Control    for down
R               to regenerate city
Y               to invert mouse Y axis

## Prerequisites & How to run:
- CMake
- C++ Compiler: MinGW-W64 or compatible
Dependencies:
- Wolf Framework
- GLEW, GLFW, GLM, stb_image, Assimp

Clone the repository, Open project in IDE.
Ctrl-Shift-P and click CMake:Configure
Choose MinGW-W64 compiler and click Configure
Then choose to run ./sample.exe file


## Requirements
For my assignment I chose to do:
- All core requirements
- Road Generation
- Lighting
- Building shapes
- Some extra stuff as well

Edits to wolf;
dfd

- Added an automatic Texture Array initialization by just adding file names [W_Texture],[W_TextureManager]
- Changed [W_App] to allow locking the cursor

## Explanation of my thought process and choices:

Attached next to this readme file is the Cityscape Documented pdf where I go more in depth about my choices and ideas throughout the whole process. Towards the end the CityScape Documented pdf ends up being more of a To-DO list rather than documenting the process. I would have liked to have more time to document the process and explain the thought process behind each choice.

Manually coded a Voronoi diagram using videos explaining it such as: 
[A Mathematical Guide to Social Distancing | Voronoi Diagrams](https://www.youtube.com/watch?v=lmbegJm4EpA&ab_channel=TwoAngles)

After getting the voronoi diagram to work, I wanted to add more procedurality so I would grab any polygon that's too big and split it into smaller ones. This is done by finding the largest edge and splitting it in half.

After getting a nice polygon, I would sweep across it cutting it into grid cells. This would result in building shapes. I used a simple grid system to do this. I also added some randomness to the grid to make it more interesting.

Then now that I have these squares, I would create the walls using them.

To add roofs to these potentially complex polygon shapes, I used a fan triangulation algorithm that would pick a vertex and then create a triangle from that vertex to the two vertices that are next in line. Mimicing a fan shape, this would create decent roof shape.

Seeing as the voronoi diagram is a set of points, it seemed like a good idea to use these points to create the roads. Because most of them are usually connecting, I would leave some padding and create a cross walk to give the roads a cohesive feel to them. And ofcourse add lines along the verticies using GL_Lines copying what a road with lines would look like.

Next up is creating UV maps for each building, this is where I faced the biggest issue with my fan based building roofs. Because the fan triangulation algorithm creates triangles that are not necessarily in order, it makes it hard to create a UV map. I ended up using a simple algorithm that would create a UV map based on the minimum x, minimum y and maximum x, maximum y of the building. This would give the texture the same orientation as the building and not be stretched. I also added some randomness to the UV map to make it more interesting.

Faced some trouble with the Wolf Texture manager, because it only takes .tga files. And most of which when cropped would lose their padding and come tilted and weird. The issue ended up getting fixed once you convert the code to take 32bit files.

I used Array Textures to store the UV maps for each building. This was a good idea because it allowed me to easily switch between different textures and also made it easy to add more textures in the future. Not to mention it made it easy to add more buildings with different textures.

I used this Array Texture to create red window variations of 3 building textures, 
so that the buildings would light up at night. And I also added normals to each building so that they would react to light.

After that I added a Day & Night cycle to the scene. This was done by changing the light color and intensity based on the time of day and have the sun rotate in the day time.

Finally I added a ground, sidewalks, building districts and special buildings to give the scene more depth and interest. This was done by creating different colourschemes and shapes for each district and giving the industrial district the special building that is stacked up into a tower shape.

I hope that the code, this and the pdf document will be enough to give you the jist of the whole project.

## More Screenshots:

![s2.png](https://github.com/Bbrizly/Cityscape/blob/main/data/screenshots/s2.png)
![s3.png](https://github.com/Bbrizly/Cityscape/blob/main/data/screenshots/s3.png)
![s5.png](https://github.com/Bbrizly/Cityscape/blob/main/data/screenshots/s5.png)
![s6.png](https://github.com/Bbrizly/Cityscape/blob/main/data/screenshots/s6.png)
![s1.png](https://github.com/Bbrizly/Cityscape/blob/main/data/screenshots/s1.png)
