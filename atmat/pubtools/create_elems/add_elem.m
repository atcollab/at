function [ ring2 ] = add_elem( ring, Elem, index )
%[ ring2 ] = add_elem( ring, Elem, index )
%add_elem adds an element to a ring at position index and outputs a new one
ringpart1 = ring(1:index);
ringpart2 = ring(index+1:end);

% ring2={ringpart1{:}; Elem; ringpart2{:}};
 ring2=[ringpart1;{Elem};ringpart2];