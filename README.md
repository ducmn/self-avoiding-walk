# self-avoiding-walk
This is my implementation of the [self-avoiding walk (SAW)](https://en.wikipedia.org/wiki/Self-avoiding_walk) in Python.

Basically, these are approaches to count how many different ways are there to "walk" a certain number of steps in a grid, without steping on a node where you have already tread. So there are 4 ways to walk one step (up, down, left, right), 12 ways to walk two steps (up left, up up, up right, right up, right right, right down, down right, down down, down left, left down, left left, left up), and so on. As you can imagine, the problem gets spicy very quickly.

This repo shows a few ways of doing this:
- Recursive: Just check every single possible directions for a valid step; very computationaly expensive
- Simple-sampling: Choose a random direction to go in every time (even if it's an invalid move), count the number of successfully completed "walks" in a number of attempts, then extrapolate to estimate how many total walks there could be; gives OK estimates, but moderately computationaly expensive
- Rosenbluth-sampling: Some complicated maths to extrapolate the number of total walks; gives acceptable estimates, but is really fast