from PIL import Image
import os

def makeGIF(name, count):
    images = []
    for i in range(0, count):
        images.append(Image.open(f"{name}{i}.png"))
    images[0].save(f'{name}.gif', save_all=True, append_images=images[1:], duration=100, loop=0)

def deletePNG(name, count):
    for i in range(0, count):
        file = f"{name}{i}.png"
        os.remove(file)

def GIF(name, count) -> None:
    makeGIF(name, count)
    deletePNG(name, count)