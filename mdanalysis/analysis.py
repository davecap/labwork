#!/usr/bin/python
# -*- coding: utf-8 -*-

class Analysis(object):
    """ Base analysis class """
    def __init__(self):
        raise NotImplementedError()
    
    def prepare(self, ref, trj):
        raise NotImplementedError()
    
    def process(self, frame):
        raise NotImplementedError()
        
    def results(self):
        raise NotImplementedError()