# Ant Simulation using Tkinter

"""
Overview:
- Uses Tkinter Canvas for graphics instead of Pygame
- Ants are represented as oval shapes on the canvas
- Animation uses the after() method for continuous updates

Ideas for Enhancement:
- Food: Add food sources that ants can collect and bring back to a nest.
- Pheromones: Implement pheromone trails that ants follow to find food and the nest.
- Different ant types: Create different types of ants with specialized roles (e.g., workers, soldiers).
- Obstacles: Add obstacles that ants must navigate around.
- Improved movement: Make the ant movement more realistic using steering behaviors or other algorithms.
"""

import tkinter as tk
import random

# Screen dimensions
screen_width = 800
screen_height = 600

# Ant class
class Ant:
    def __init__(self, canvas, x, y):
        self.canvas = canvas
        self.x = x
        self.y = y
        self.size = 5
        # Create oval on canvas
        self.id = canvas.create_oval(
            x - self.size, y - self.size,
            x + self.size, y + self.size,
            fill='black'
        )
        self.direction = random.randint(0, 359)

    def update(self):
        # Move the ant
        dx = 2 * random.randint(-1, 1)
        dy = 2 * random.randint(-1, 1)
        
        # Update position with boundary checking
        self.x = max(self.size, min(self.x + dx, screen_width - self.size))
        self.y = max(self.size, min(self.y + dy, screen_height - self.size))
        
        # Update position on canvas
        self.canvas.coords(
            self.id,
            self.x - self.size, self.y - self.size,
            self.x + self.size, self.y + self.size
        )

class AntSimulation:
    def __init__(self, root):
        self.root = root
        self.root.title("Ant Simulation")
        
        # Create canvas
        self.canvas = tk.Canvas(
            root, 
            width=screen_width, 
            height=screen_height, 
            bg='white'
        )
        self.canvas.pack()
        
        # Create ants
        self.ants = []
        for i in range(10):
            x = random.randint(0, screen_width)
            y = random.randint(0, screen_height)
            self.ants.append(Ant(self.canvas, x, y))
        
        # Start animation
        self.running = True
        self.animate()
        
        # Bind close event
        self.root.protocol("WM_DELETE_WINDOW", self.on_close)
    
    def animate(self):
        if self.running:
            # Update all ants
            for ant in self.ants:
                ant.update()
            
            # Schedule next frame (approximately 60 FPS)
            self.root.after(16, self.animate)
    
    def on_close(self):
        self.running = False
        self.root.destroy()

def main():
    root = tk.Tk()
    app = AntSimulation(root)
    root.mainloop()

if __name__ == "__main__":
    main()
