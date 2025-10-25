# Ant Simulation using Matplotlib
# Install: pip install matplotlib

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random

# Screen dimensions
screen_width = 800
screen_height = 600

# Ant class
class Ant:
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.direction = random.randint(0, 359)

    def update(self):
        # Move the ant
        dx = 2 * random.randint(-1, 1)
        dy = 2 * random.randint(-1, 1)
        
        # Update position with boundary checking
        self.x = max(0, min(self.x + dx, screen_width))
        self.y = max(0, min(self.y + dy, screen_height))

def main():
    # Create figure and axis
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_xlim(0, screen_width)
    ax.set_ylim(0, screen_height)
    ax.set_title("Ant Simulation")
    ax.set_facecolor('white')
    
    # Create ants
    ants = []
    for i in range(10):
        x = random.randint(0, screen_width)
        y = random.randint(0, screen_height)
        ants.append(Ant(x, y))
    
    # Create scatter plot for ants
    scatter = ax.scatter([], [], c='black', s=50)
    
    def update_frame(frame):
        # Update all ants
        for ant in ants:
            ant.update()
        
        # Update scatter plot
        positions = [(ant.x, ant.y) for ant in ants]
        scatter.set_offsets(positions)
        return scatter,
    
    # Create animation
    ani = animation.FuncAnimation(
        fig, update_frame, interval=16, blit=True, cache_frame_data=False
    )
    
    plt.show()

if __name__ == "__main__":
    main()
