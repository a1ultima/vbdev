#! /usr/bin/python3

class Tree:
    def __init__ (self, weight, children):
        self.weight = weight
        self.children = children [:]

    def distances (self, curDistance = .0, acc = None):
        if acc is None: acc = []
        for child in self.children:
            child.distances (self.weight + curDistance, acc)
        return acc

    def collapse (self, limit):
        self.children = [child.collapse (limit) for child in self.children]
        distances = self.distances (-self.weight)
        avg = sum (distances) / len (distances)
        if avg > limit: return self
        return Node (self.weight, ''.join (self.descendants () ) )

    def descendants (self):
        descendants = []
        for child in self.children:
            descendants.extend (child.descendants () )
        return descendants

    def __repr__ (self):
        return '({}):{}'.format (','.join (str (child) for child in self.children), self.weight)

class Node:
    def __init__ (self, weight, name):
        self.weight = weight
        self.name = name

    def distances (self, curDistance, acc):
        acc.append (curDistance + self.weight)

    def collapse (self, limit):
        return self

    def descendants (self):
        return [self.name]

    def __repr__ (self):
        return '{}:{}'.format (self.name, self.weight)

class Stack (list):
    def pop (self):
        e = self [0]
        del self [0]
        return e

    def push (self, e):
        self.insert (0, e)

def parse (tree):
    buff = ''
    stack = Stack ()
    while True:
        c = tree [0]
        if c == ';': break
        tree = tree [1:]

        if c == '(':
            stack.push (c)
            continue

        if c in ':,':
            if buff: stack.push (buff)
            buff = ''
            continue

        if c == ')':
            if buff: stack.push (buff)
            buff = ''
            popped = ''
            children = []
            while True:
                weight = stack.pop ()
                if weight == '(': break
                weight = float (weight)
                child = stack.pop ()
                if isinstance (child, Tree):
                    child.weight = weight
                else:
                    child = Node (weight, child)
                children.append (child)
            stack.push (Tree (0, children) )
            continue

        buff += c

    return stack.pop ()

t = parse ('((A:0.9,(B:0.2,C:0.3):0.3,(E:0.05,F:0.08):0.1):0.6);')
print ('Input tree is {}'.format (t) )
for limit in range (1, 6):
    limit = limit / 10
    print ('Collapsed by {} is {}'.format (limit, t.collapse (limit) ) )