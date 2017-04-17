/*
 * Copyright (C) 2012 Jordan Fish <fishjord at msu.edu>
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

package edu.msu.cme.rdp.graph.hash;

/**
 *
 * @author wangqion
 */
public interface Hash {

    public long eatRight(long hashvalue, int c);
    public long eatLeft(long hashvalue, int c);

    //public long getInitialHashvalue(String s);
    public long updateLeft(long hashvalue , int outint, int inint );
    public long updateRight(long hashvalue , int outint, int inint );

    public long replaceRight(long hashvalue, int oldC, int newC) ;
    public long replaceLeft( long hashvalue, int oldC, int newC);
}
