/*
 * Copyright (C) 2015 leotift
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.msu.cme.rdp.primerdesign.deadcode;

/**
 *
 * @author leotift
 */
public class RDPStringBuilder implements Cloneable{
    
    protected StringBuilder stringBuilder = new StringBuilder();
          
    protected RDPStringBuilder() {}
       
    public StringBuilder getStringBuilder() {
        return this.stringBuilder;   
    }
    
    @Override
    public Object clone() {
      Object clone = null;
      
      try {
         clone = super.clone();
         
      } catch (CloneNotSupportedException e) {
          System.out.println("Cloning sequence problem");
      }
      
      return clone;
    }
  
}
