#Ref from generative AI

import pygame
import random
import math

pygame.init()

SCREEN_WIDTH = 800
SCREEN_HEIGHT = 600
screen = pygame.display.set_mode((SCREEN_WIDTH, SCREEN_HEIGHT))
pygame.display.set_caption("Ant Simulation with Food")

BLACK = (0, 0, 0)
WHITE = (255, 255, 255)
RED = (255, 0, 0)
GREEN = (0, 200, 0)
BLUE = (0, 0, 255)

class Nest(pygame.sprite.Sprite):
    def __init__(self, x, y):
        super().__init__()
        self.image = pygame.Surface((24, 24))
        self.image.fill(RED)
        self.rect = self.image.get_rect(center=(x, y))
        self.food_collected = 0


class Food(pygame.sprite.Sprite):
    def __init__(self, x, y):
        super().__init__()
        self.image = pygame.Surface((8, 8))
        self.image.fill(GREEN)
        self.rect = self.image.get_rect(center=(x, y))
        self.collected = False


class Ant(pygame.sprite.Sprite):
    def __init__(self, x, y, nest):
        super().__init__()
        self.image = pygame.Surface((6, 6))
        self.image.fill(BLACK)
        self.rect = self.image.get_rect(center=(x, y))

        self.pos = pygame.math.Vector2(self.rect.center)
        self.speed = 1.8
        self.state = "searching"
        self.target_food = None
        self.nest = nest

    def update(self, food_group):
        if self.state == "searching":
            jitter = pygame.math.Vector2(random.uniform(-1.8, 1.8), random.uniform(-1.8, 1.8))
            self.pos += jitter

            self.pos.x = max(0, min(self.pos.x, SCREEN_WIDTH))
            self.pos.y = max(0, min(self.pos.y, SCREEN_HEIGHT))
            self.rect.center = (round(self.pos.x), round(self.pos.y))

            for food in list(food_group):
                if self.rect.colliderect(food.rect):
                    self.state = "returning"
                    self.target_food = food
                    food.collected = True
                    food_group.remove(food)
                    break

        elif self.state == "returning":
            nest_vec = pygame.math.Vector2(self.nest.rect.center)
            direction = (nest_vec - self.pos)
            distance = direction.length()
            if distance > 0:
                direction = direction.normalize()
                self.pos += direction * self.speed

            self.rect.center = (round(self.pos.x), round(self.pos.y))

            if self.rect.colliderect(self.nest.rect):
                self.nest.food_collected += 1
                self.state = "searching"
                self.target_food = None
                away = pygame.math.Vector2(random.uniform(-5, 5), random.uniform(-5, 5))
                self.pos += away
                self.rect.center = (round(self.pos.x), round(self.pos.y))

def main():
    clock = pygame.time.Clock()

    nest = Nest(SCREEN_WIDTH // 2, SCREEN_HEIGHT // 2)

    ants = pygame.sprite.Group()
    for _ in range(20):
        ants.add(Ant(nest.rect.centerx, nest.rect.centery, nest))

    food_group = pygame.sprite.Group()
    for _ in range(30):
        angle = random.uniform(0, 2 * math.pi)
        radius = random.uniform(30, 120)
        x = int(nest.rect.centerx + radius * math.cos(angle))
        y = int(nest.rect.centery + radius * math.sin(angle))

        x = max(10, min(SCREEN_WIDTH - 10, x))
        y = max(10, min(SCREEN_HEIGHT - 10, y))

        food_group.add(Food(x, y))

        font = pygame.font.SysFont(None, 20)

    running = True
    while running:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False

            if event.type == pygame.MOUSEBUTTONDOWN and event.button == 1:
                mx, my = event.pos
                food_group.add(Food(mx, my))

        ants.update(food_group)

        screen.fill(WHITE)

        screen.blit(nest.image, nest.rect)

        food_group.draw(screen)
        ants.draw(screen)

        remaining = len(food_group)
        hud = font.render(f"Food remaining: {remaining}    Food at nest: {nest.food_collected}", True, BLUE)
        screen.blit(hud, (10, 10))

        pygame.display.flip()
        clock.tick(60)

    pygame.quit()


if __name__ == "__main__":
    main()
