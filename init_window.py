#!/opt/local/bin/python
from Tkinter import *  
import ttk as ttk
import tkMessageBox
import numpy as np
from weight_inverse_distance import *

M = np.matrix('-1. 3. -3. 1.; 3. -6. 3. 0.; -3. 3. 0. 0.; 1. 0. 0. 0.')

gOnce = False

class Window:
	canvas = None
	mode = None #indicator of one of the three modes.
	root = None
	selected = None
	popup = None
	trace = None
	transforms = {}
	
	def __init__(self, parent):
		self.root = parent
		mainWindow = ttk.Frame(parent)	
		self.init_UI(mainWindow)	
		mainWindow.grid()
		return
		
	def init_UI(self, parent):
		menubar = ttk.Frame(parent, relief=GROOVE, borderwidth=1, width=600, height=28)
		self.init_menubar(menubar)
		
		self.canvas = Canvas(parent, width=600, height=400, bd=2, cursor='pencil', relief=SUNKEN)
		self.canvas.bind("<Button-1>", self.onclick_handler)
		self.canvas.bind("<Shift-B1-Motion>", self.on_shift_mouse_handler)
		self.canvas.bind("<Control-B1-Motion>", self.on_control_mouse_handler)
#		self.canvas.bind("<Double-Button-1>", self.on_double_click_handler)
		self.canvas.bind("<B1-Motion>", self.on_motion_handler)
		self.canvas.bind("<ButtonRelease-1>", self.onrelease_handler)
		self.canvas.grid()
		return	
		
	def onclick_handler(self, event):
		mode = self.mode.get()
		if mode == 0:
			controls = self.canvas.find_withtag('controls')
			if len( controls ) == 4:
				self.canvas.delete( controls[0] )
			
			r= 3
			x0, x1, y0, y1 = event.x-r, event.x+r, event.y-r, event.y+r
			self.canvas.create_oval(x0, y0, x1, y1, fill='blue', outline='blue', tags='controls')
			
			if len( self.canvas.find_withtag('controls') ) == 4:
				self.redraw_bezier_curve( )
				if len(self.canvas.find_withtag('handles')) > 0:
					self.redraw_handle_affected_curve()
					self.redraw_approximated_bezier()

		elif mode == 1:
			r = 3
			x0, x1, y0, y1 = event.x-r, event.x+r, event.y-r, event.y+r
			handle = self.canvas.create_rectangle(x0, y0, x1, y1, fill='red', outline='red', tags='handles')
			self.transforms[handle] = np.identity(3).reshape(-1)
			self.popup_handle_editor( handle )
			if len( self.canvas.find_withtag('controls') ) == 4:
				self.redraw_handle_affected_curve()
				self.redraw_approximated_bezier()

		elif mode == 2:
			overlaps = self.canvas.find_overlapping(event.x-3, event.y-3, event.x+3, event.y+3)
			sels = set(overlaps) & set(self.canvas.find_withtag('handles')) 
			if len(sels) > 0:
				self.selected = sels.pop()
				coord = self.canvas.coords( self.selected )
				self.trace = [(coord[0]+coord[2])/2, (coord[1]+coord[3])/2]
				self.popup_handle_editor(self.selected)
				
#	def on_double_click_handler(self, event):		
#		if self.mode.get() == 2:
#			overlaps = self.canvas.find_overlapping(event.x-3, event.y-3, event.x+3, event.y+3)
#			sels = set(overlaps) & set(self.canvas.find_withtag('handles')) 
#			if len(sels) > 0:
#				self.popup_handle_editor(sels.pop())
	
	def on_motion_handler(self, event):
		if self.selected != None and self.popup != None: 
			h = self.selected
			trans_x, trans_y = event.x-self.trace[0], event.y-self.trace[1]
			self.trace = [event.x, event.y]
			
			window_name = self.canvas.itemcget( self.popup, 'window' )
			p = self.root.nametowidget( window_name )

			labelFrame = p.winfo_children()[0]
			entry_transX = labelFrame.winfo_children()[2]
			entry_transY = labelFrame.winfo_children()[5]
			
			ori_x, ori_y = self.transforms[h][2], self.transforms[h][5]
			new_x, new_y = trans_x+ori_x, trans_y+ori_y
			entry_transX.delete(0,END)
			entry_transX.insert(0,new_x)
			entry_transY.delete(0,END)
			entry_transY.insert(0,new_y)
			self.transforms[h][2], self.transforms[h][5] = new_x, new_y
			
#			self.canvas.move(self.selected, event.x-ori_x, event.y-ori_y)
			if len( self.canvas.find_withtag('controls') ) == 4:
				self.redraw_handle_affected_curve()
				self.redraw_approximated_bezier()				
		return
		
	def on_shift_mouse_handler(self, event):
		if self.selected != None and self.popup != None: 
			print 'hello world'
		
	def on_control_mouse_handler(self, event):
		if self.selected != None and self.popup != None: 
			print 'control'
	
	def redraw_bezier_curve(self):
		self.canvas.delete( 'original_bezier' )
		cps = self.canvas.find_withtag('controls')
		cps = [self.canvas.bbox(x) for x in cps]
		cps = [((x[0]+x[2])/2, (x[1]+x[3])/2) for x in cps]
		cps = np.mat(cps)
		known = np.asarray( M*cps ).reshape(4, -1)
		
		ps = [] 
		for t in range(0, 101):
			p = np.dot( np.asarray( [(float(t)/100)**3, (float(t)/100)**2, float(t)/100, 1] ), known )
			ps = ps + [p[0], p[1]]	
		self.canvas.create_line(ps, smooth=True, width=2, tags='original_bezier')
		
	def redraw_handle_affected_curve(self):
		self.canvas.delete( 'handle_affected' )
		pts = self.canvas.find_withtag( 'original_bezier' )
		pts = self.canvas.coords( pts )
		
		tps = []
		for i in range(0, len(pts)/2):
			p = np.asarray([pts[i*2], pts[i*2+1], 1.])
			m = np.zeros(9)
			w = {}
			for h in self.canvas.find_withtag('handles'):
			# coefficient decided by square of distance
				pos = self.canvas.coords(h)
				w[h] = 1 / ( (pos[0]/2+pos[2]/2-p[0])**2 + (pos[1]/2+pos[3]/2-p[1])**2 )
				
			for h in self.canvas.find_withtag('handles'):
				m = m + self.transforms[h]*( w[h] / sum(w.values()) )
			
			p = np.dot( m.reshape(3, 3), p.reshape(3,-1) ).reshape(-1)
			tps = tps + [p[0], p[1]]

		self.canvas.create_line(tps, smooth=True, width=2, fill='magenta', tags='handle_affected')
		
	def redraw_approximated_bezier(self):
		self.canvas.delete( 'approximated' )
		self.canvas.delete( 'new_controls' )
		
		handles = []
		H_set = self.canvas.find_withtag('handles')
		for i in range( len( H_set ) ):
			pos = self.canvas.coords( H_set[i] )
			pos = [(pos[0]+pos[2])/2, (pos[1]+pos[3])/2, 1.]
			handles.append(pos)
		
		cps = self.canvas.find_withtag('controls')
		cps = [self.canvas.bbox(x) for x in cps]
		cps = [((x[0]+x[2])/2, (x[1]+x[3])/2, 1) for x in cps]
		cps = np.mat(cps)
		
		P_prime = np.zeros(12).reshape(3,4)
		'''
		global gOnce
		
		import csv
		outname = 'integration_accuracy.csv'
		def uniquepath( path ):
			import os
			i = 1
			result = path
			while os.path.exists( result ):
				split = os.path.splitext( path )
				result = split[0] + ( ' %d' % i ) + split[1]
				i += 1
			return result
		
		writer = csv.writer( open( uniquepath( 'integration_accuracy.csv' ), 'w' ) )
		'''
		for i in range( len( H_set ) ):

			T_i = np.mat( np.asarray(self.transforms[H_set[i]]).reshape(3,3) )
			W_i = precompute_W_i( handles, i, cps, M, 0., 1., 50 )
			
			#for num_samples in xrange(1,501):
			#	writer.writerow( [ num_samples, i ] + list( precompute_W_i( handles, i, cps, M, 0., 1., num_samples ).flat ) )
			
			inv_A = precompute_inv_A( M, 0., 1., 50 )
			P_prime = P_prime + T_i * (cps.T) * M * np.mat(W_i) * inv_A
		
		#del writer
		
		for pp in P_prime.T:
			pp = np.asarray(pp).reshape(3)
			r= 3
			x0, x1, y0, y1 = pp[0]-r, pp[0]+r, pp[1]-r, pp[1]+r
			self.canvas.create_oval(x0, y0, x1, y1, fill='cyan', outline='cyan', tags='new_controls')

		
		known = np.asarray( M * P_prime.T ).reshape(4, -1)	
		ps = [] 
		for t in range(0, 101):
			p = np.dot( np.asarray( [(float(t)/100)**3, (float(t)/100)**2, float(t)/100, 1] ), known )
			ps = ps + [p[0], p[1]]	
		self.canvas.create_line(ps, smooth=True, width=2, fill='green', tags='approximated')	 

	
	def onrelease_handler(self, event):
		self.trace = []
		if self.selected == None:
			self.canvas.delete('popup') 
		else:	
			self.selected = None

	def change_mode(self):
		mode = self.mode.get()
		if mode == 0:
			self.canvas.config(cursor = 'pencil')
		elif mode == 1: 
			self.canvas.config(cursor = 'target')
		elif mode == 2:
			self.canvas.config(cursor = 'hand1')

	def del_controls(self):
		self.canvas.delete('controls')
		self.canvas.delete('original_bezier')
		self.canvas.delete('handle_affected')
		self.canvas.delete('approximated')
		
	
	def del_handles(self):
		self.canvas.delete('handles')
		self.canvas.delete('handle_affected')
		self.canvas.delete('approximated')
		self.transforms.clear()

	def clear_canvas(self):
		self.canvas.delete('all')
		self.transforms.clear()
	
	def save_transforms(self, id, vals, popup):
		
		for i in range(0, 9):
			self.transforms[id][i] = vals[i].get()
		self.remove_popup(popup)
		
		if len( self.canvas.find_withtag('controls') ) == 4:
			self.redraw_handle_affected_curve()
			self.redraw_approximated_bezier()	
		return
		
	def delete_handle(self, id, popup):
		self.canvas.delete(id)
		self.remove_popup(popup)
		
		if len( self.canvas.find_withtag('controls') ) == 4:
			self.redraw_handle_affected_curve()
			self.redraw_approximated_bezier()	
		return	
		
	def remove_popup(self, popup):
		self.canvas.delete(popup)
		
	def init_menubar(self, menubar):
		menubar.grid_propagate(0)
		menubar.grid()
		
		#A menu in Tk is a combination of a Menubutton (the title of the
		#menu) and the Menu (what drops down when the Menubutton is pressed)
			
		mb_file = ttk.Menubutton(menubar, text='file')
		mb_file.grid(column=0, row=0)
		mb_file.menu = Menu(mb_file)
		
		#Once we've specified the menubutton and the menu, we can add
		#different commands to the menu
		
		mb_file.menu.add_command(label='open')
		mb_file.menu.add_command(label='close')
		
		mb_edit = ttk.Menubutton(menubar, text='edit')
		mb_edit.grid(column=1, row=0)
		mb_edit.menu = Menu(mb_edit)
		
		self.mode = IntVar()
		mb_edit.menu.add_radiobutton(label='add control point', variable=self.mode, value=0, command=self.change_mode)
		mb_edit.menu.add_radiobutton(label='add handle', variable=self.mode, value=1, command=self.change_mode)
		mb_edit.menu.add_radiobutton(label='edit handle', variable=self.mode, value=2, command=self.change_mode)
		mb_edit.menu.add_separator()
		mb_edit.menu.add_command(label='delete controls', command=self.del_controls)
		mb_edit.menu.add_command(label='delete handles', command=self.del_handles)
		mb_edit.menu.add_command(label='clear canvas', command=self.clear_canvas)
		
		
		mb_help = ttk.Menubutton(menubar,text='help')
		mb_help.grid(column=2, row=0, padx=340, sticky=E)
		
		mb_file['menu'] = mb_file.menu
		mb_edit['menu'] = mb_edit.menu
		return 
		
	def popup_handle_editor(self, handle_id):
	
		self.canvas.delete('popup') 
		self.selected = handle_id
		
		coord = self.canvas.bbox(handle_id)
		data = self.transforms[handle_id]
		values = []
		for i in range(0, 9):
			values.append(StringVar())
			values[i].set(data[i])
		
		frame = Frame(self.root, borderwidth=1, relief=RIDGE)
		labelFrame = LabelFrame(frame, text='Transform', relief=GROOVE, borderwidth=1)
		labelFrame.grid()
		w11 = Entry(labelFrame, cursor='xterm', highlightcolor='cyan', width=6, textvariable=values[0])
		w11.grid(row=0, column=0)
		w12 = Entry(labelFrame, cursor='xterm', highlightcolor='cyan', width=6, textvariable=values[1])
		w12.grid(row=0, column=1)
		w13 = Entry(labelFrame, cursor='xterm', highlightcolor='cyan', width=6, textvariable=values[2])
		w13.grid(row=0, column=2)
		w21 = Entry(labelFrame, cursor='xterm', highlightcolor='cyan', width=6, textvariable=values[3])
		w21.grid(row=1, column=0)
		w22 = Entry(labelFrame, cursor='xterm', highlightcolor='cyan', width=6, textvariable=values[4])
		w22.grid(row=1, column=1)
		w23 = Entry(labelFrame, cursor='xterm', highlightcolor='cyan', width=6, textvariable=values[5])
		w23.grid(row=1, column=2)
		w31 = Entry(labelFrame, cursor='xterm', highlightcolor='cyan', width=6, textvariable=values[6])
		w31.grid(row=2, column=0)
		w32 = Entry(labelFrame, cursor='xterm', highlightcolor='cyan', width=6, textvariable=values[7])
		w32.grid(row=2, column=1)
		w33 = Entry(labelFrame, cursor='xterm', highlightcolor='cyan', width=6, textvariable=values[8])
		w33.grid(row=2, column=2)
		
		popup = self.canvas.create_window((coord[0]+coord[2])/2, (coord[1]+coord[3])/2, anchor=NW, width=204, height=140, window=frame, tags='popup')	
		
		btn1 = Button(frame, text='Save', command=lambda id=handle_id, vals=values, popup=popup : self.save_transforms(id, vals, popup), width=4)
		btn1.grid(row=3, column=0, sticky=SW)
		btn2 = Button(frame, text='Delete', command=lambda id=handle_id, popup=popup : self.delete_handle(id, popup), width=4)
		btn2.grid(row=3, column=0, sticky=S)
		btn3 = Button(frame, text='Cancel', command=lambda popup=popup : self.remove_popup(popup), width=4)
		btn3.grid(row=3, column=0, sticky=SE)
	
		self.popup = popup
#		return popup

	def onExit(self):
		self.quit()


def main():
  
	root = Tk()
#	root.geometry("600x400+20+20")
	app = Window(root)
	root.mainloop()	 

if __name__ == '__main__': main()
